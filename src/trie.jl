abstract AbstractTrieNode

immutable NullTrieNode <: AbstractTrieNode
end

type TrieNode{A<:Alphabet,K} <: AbstractTrieNode
   next::Vector{AbstractTrieNode}
   offsets::Vector{Vector{UInt32}}
   metadata::Vector{Vector{K}}

   function TrieNode()
      next = Vector{AbstractTrieNode}(length(alphabet(A)))
      fill!(next, NullTrieNode())
      return new(next, Vector{Vector{UInt32}}(length(alphabet(A))), 
                       Vector{Vector{K}}(length(alphabet(A))))
   end
end

type RNATrie{A<:Alphabet,K}
   root::AbstractTrieNode
   range::UnitRange

   RNATrie(range::UnitRange) = new(NullTrieNode(), range)
end

function Base.push!{A,K}( trie::RNATrie{A,K}, seq::Bio.Seq.Sequence, metadata::K )
   if isa( trie.root, NullTrieNode )
      trie.root = TrieNode{A,K}()
   end
   for i in 1:length(seq)-first(trie.range)+1
      const range = i+last(trie.range)-1 > length(seq) ? (first(trie.range):length(seq)-i+1) : trie.range
      push!( trie.root, seq, trie.range, i, metadata )
   end
end

function Base.push!{A,K}( node::TrieNode{A,K}, seq::Bio.Seq.Sequence, 
                          range::UnitRange, idx::Int, metadata::K; 
                          curdepth::Int=1 )
   if curdepth > last(range) || idx > length(seq) || isambiguous( seq[idx] )
      return
   end
   nucidx = Bio.Seq.encode( A, seq[idx] ) + 1
   if isa( node.next[nucidx], NullTrieNode )
      node.next[nucidx] = TrieNode{A,K}()
   end
   if curdepth in range
      if !isdefined( node.offsets, nucidx )
         node.offsets[nucidx]  = Vector{UInt32}()
         node.metadata[nucidx] = Vector{K}()
      end
      push!( node.offsets[nucidx],  convert(UInt32, idx) )
      push!( node.metadata[nucidx], metadata )
   end
   push!( node.next[nucidx], seq, range, idx + 1, metadata, curdepth=curdepth + 1 )
end

nodecount{A,K}( trie::RNATrie{A,K} ) = nodecount( trie.root ) - 1
nodecount( node::NullTrieNode ) = 0

function nodecount{A,K}( node::TrieNode{A,K} )
   cnt = 1
   for n in node.next
      cnt += nodecount( n )
   end
   cnt
end

immutable DuplexTrieSeq
   name::String
   seq::Bio.Seq.Sequence
   strand::Bool
   genomeoffset::Int
   length::Int
end

type DuplexTrie{A <: Alphabet,K <: Integer}
   fwd::RNATrie{A,K}
   rev::RNATrie{A,K}
   entry::Vector{DuplexTrieSeq}
   range::UnitRange
   
   function DuplexTrie( seq::Bio.Seq.Sequence, range::UnitRange; 
                        seqname::String="chr", genomepos::Int=1, strand::Bool=true )
      fwd = RNATrie{A,K}( range )
      rev = RNATrie{A,K}( range )
      entry = DuplexTrieSeq( seqname, seq, strand, genomepos-1, length(seq) )
      push!( fwd, seq, one(K) )
      push!( rev, reverse(seq), one(K) )
      return new( fwd, rev, DuplexTrieSeq[entry], range )
   end

   function DuplexTrie( left::Bio.Seq.Sequence, right::Bio.Seq.Sequence, range::UnitRange; 
                        left_seqname::String="chr", right_seqname::String="chr",
                        left_genomepos::Int=1, right_genomepos::Int=1,
                        left_strand::Bool=true, right_strand::Bool=true )
      fwd = RNATrie{A,K}( range )
      rev = RNATrie{A,K}( range ) 
      left_entry  = DuplexTrieSeq( left_seqname, left, left_strand, left_genomepos-1, length(left) )
      right_entry = DuplexTrieSeq( right_seqname, right, right_strand, right_genomepos-1, length(right) )
      push!( fwd, left, one(K) )
      push!( rev, reverse(right), convert(K, 2) )
      return new( fwd, rev, DuplexTrieSeq[left_entry, right_entry], range )
   end

   function DuplexTrie( seqs::Vector{Bio.Seq.Sequence}, range::UnitRange;
                        seqnames::Vector{String}=String["chr" for i in 1:length(seqs)],
                        genomepos::Vector{Int}=ones(Int, length(seqs)),
                        strands::Vector{Bool}=trues(length(seqs)))
      fwd = RNATrie{A,K}( range )
      rev = RNATrie{A,K}( range )
      genomepos -= 1
      entries = DuplexTrieSeq[ DuplexTrieSeq( seqnames[i], seqs[i], strands[i], genomepos[i], length(seqs[i]) ) for i in 1:length(seqs) ]
      @assert typemax(K) >= length(seqs)
      for i in 1:length(seqs)
         push!( fwd, seqs[i], convert(K, i) )
         push!( rev, reverse(seqs[i]), convert(K, i) )
      end
      return new( fwd, rev, entries, range )
   end
end

function traverse{A,K}( trie::DuplexTrie{A,K}, foldrange::UnitRange; 
                        bulge_max::Int=zero(Int), mismatch_max::Int=zero(Int),
                        raw_output::String="")

   const pairs_idx      = index(A, pairs)
   const bulges_idx     = index(A, gaps)
   const mismatches_idx = index(A, mismatches)

   const duplex         = RNADuplex()
   const deprange       = trie.range
   const energy_max     = 0.01

   const intervals      = DuplexCollection{K}()

   if raw_output != ""
      const io = open( raw_output, "w" )
      const stream = ZlibDeflateOutputStream( io )
      const output = true
   else
      const output = false
   end

   @inline function traverse{A,K}( fwd::TrieNode{A,K}, rev::TrieNode{A,K},
                                   fdepth::Int64, rdepth::Int64,
                                   bulge_n::Int64, mismatch_n::Int64,
                                   from_bulge::Bool, bulge_left::Bool )

      # recurse through valid watson-crick pairs
      for (l,r) in pairs_idx
         if isa( fwd.next[l], TrieNode{A,K} ) && isa( rev.next[r], TrieNode{A,K} )
            push!( duplex, convert(RNAPair, onehot(l), onehot(r)) )
            if output
               tab_write( stream, (depth, energy(duplex)))
            end
            if fdepth in deprange && rdepth in deprange && energy(duplex) < trie.range.start*-1
               for (ix,i) in enumerate(fwd.offsets[l]), (jx,j) in enumerate(rev.offsets[r])
                     const k = revoffset( j, trie.entry[ rev.metadata[r][jx] ].length )
                     fwd.metadata[l][ix] == rev.metadata[r][jx] && !((k - i) + 1 in foldrange) && continue
                     const newdup = deepcopy(duplex)
                     const fwd_entry = trie.entry[ fwd.metadata[l][ix] ]
                     const rev_entry = trie.entry[ rev.metadata[r][jx] ]
                     const fg = fwd_entry.genomeoffset
                     const rg = rev_entry.genomeoffset
                     push!( intervals, DuplexInterval( Interval(fwd_entry.name, fg+(i-fdepth+1), fg+i, '?', fwd.metadata[l][ix]), 
                                                       Interval(rev_entry.name, rg+k, rg+(k+rdepth-1), '?', rev.metadata[r][jx]),
                                                       newdup ) )
               end
            end
            traverse( fwd.next[l], rev.next[r],
                       fdepth + 1, rdepth + 1,
                       bulge_n, mismatch_n,
                       false, false )
         end
      end
      # recurse through bulges
      if bulge_n < bulge_max && fdepth > 1 && rdepth > 1 && duplex.energy[end] < energy_max
         for (l,r) in bulges_idx
            if l == 0
               (from_bulge && !bulge_left) && continue # only bulge one way
               if isa( rev.next[r], TrieNode{A,K} )
                  push!( duplex, convert(RNABulge, zero(UInt8), onehot(r)) )
                  traverse( fwd, rev.next[r],
                             fdepth, rdepth + 1,
                             bulge_n + 1, mismatch_n,
                             true, true )
               end
            elseif r == 0
               (from_bulge && bulge_left) && continue 
               if isa( fwd.next[l], TrieNode{A,K} )
                  push!( duplex, convert(RNABulge, onehot(l), zero(UInt8)) )
                  traverse( fwd.next[l], rev,
                             fdepth + 1, rdepth,
                             bulge_n + 1, mismatch_n,
                             true, false )
               end
            end
         end
      end
      # recurse through mismatches
      if mismatch_n < mismatch_max && fdepth > 1 && rdepth > 1 && 
         !from_bulge && duplex.energy[end] < energy_max
         for (l,r) in mismatches_idx
            if isa( fwd.next[l], TrieNode{A,K} ) && isa( rev.next[r], TrieNode{A,K} )
             push!( duplex, convert(RNAMismatch, onehot(l), onehot(r)) )
             traverse( fwd.next[l], rev.next[r],
                        fdepth + 1, rdepth + 1,
                        bulge_n, mismatch_n + 1,
                        from_bulge, bulge_left )
            end
         end
      end

      pop!( duplex ) # clean up this level
      return
   end

   traverse( trie.fwd.root, trie.rev.root, 
             one(Int), one(Int),
             zero(Int), zero(Int),
             false, false )

   if output
      close(stream)
      close(io)
   end

   intervals
end


abstract AbstractTrieNode

typealias Alphabet Bio.Seq.Alphabet

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
   if length(seq) <= last(trie.range)
      push!( trie.root, seq, trie.range, 1, metadata )
   else
      for i in 1:length(seq)-last(trie.range)+1
         push!( trie.root, seq, trie.range, i, metadata )
      end
   end
end

function Base.push!{A,K}( node::TrieNode{A,K}, seq::Bio.Seq.Sequence, 
                          range::UnitRange, idx::Int, metadata::K; 
                          curdepth::Int=1 )
   if curdepth > last(range) || isambiguous( seq[idx] )
      return
   end
   nucidx = Bio.Seq.encode( A, seq[idx] ) + 1
   if isa( node.next[nucidx], NullTrieNode )
      node.next[nucidx] = TrieNode{A,K}()
   end
   if curdepth in range
      #=if length(node.offsets) == 0
         node.offsets  = Vector{Vector{UInt32}}(length(alphabet(A)))
         node.metadata = Vector{Vector{K}}(length(alphabet(A)))
      end=#
      if !isdefined( node.offsets, nucidx )
         node.offsets[nucidx]  = Vector{UInt32}()
         node.metadata[nucidx] = Vector{K}()
      end
      push!( node.offsets[nucidx],  convert(UInt32, idx) )
      push!( node.metadata[nucidx], metadata )
   end
   push!( node.next[nucidx], seq, range, idx + 1, metadata, curdepth=curdepth + 1 )
end


revoffset( x::Int, seq::BioSequence ) = revoffset( x, length(seq) )
revoffset( x::Int, len::Int ) = len - x + 1

const DNAPAIRS = [(DNA_A, DNA_T),   (DNA_T, DNA_A),
                  (DNA_G, DNA_C),   (DNA_C, DNA_G),
                  (DNA_G, DNA_T),   (DNA_T, DNA_G)]

const DNAGAPS  = [(DNA_A, DNA_Gap), (DNA_Gap, DNA_A),
                  (DNA_C, DNA_Gap), (DNA_Gap, DNA_C),
                  (DNA_G, DNA_Gap), (DNA_Gap, DNA_G),
                  (DNA_T, DNA_Gap), (DNA_Gap, DNA_T)]

const DNAMISMATCH = [(DNA_A, DNA_A), (DNA_A, DNA_C),
                     (DNA_A, DNA_G), (DNA_C, DNA_A),
                     (DNA_C, DNA_C), (DNA_C, DNA_T),
                     (DNA_G, DNA_A), (DNA_G, DNA_G),
                     (DNA_T, DNA_C), (DNA_T, DNA_T)]

const RNAPAIRS = [(RNA_A, RNA_U), (RNA_U, RNA_A),
                  (RNA_G, RNA_C), (RNA_C, RNA_G),
                  (RNA_G, RNA_U), (RNA_U, RNA_G)]

const RNAGAPS  = [(RNA_A, RNA_Gap), (RNA_Gap, RNA_A),
                  (RNA_C, RNA_Gap), (RNA_Gap, RNA_C),
                  (RNA_G, RNA_Gap), (RNA_Gap, RNA_G),
                  (RNA_U, RNA_Gap), (RNA_Gap, RNA_U)]

const RNAMISMATCH = [(RNA_A, RNA_A), (RNA_A, RNA_C),
                     (RNA_A, RNA_G), (RNA_C, RNA_A),
                     (RNA_C, RNA_C), (RNA_C, RNA_U),
                     (RNA_G, RNA_A), (RNA_G, RNA_G),
                     (RNA_U, RNA_C), (RNA_U, RNA_U)]

typealias PairsType Vector{Tuple{Bio.Seq.Nucleotide, Bio.Seq.Nucleotide}}

pairs{n}(::Type{Bio.Seq.DNAAlphabet{n}})      = DNAPAIRS
pairs{n}(::Type{Bio.Seq.RNAAlphabet{n}})      = RNAPAIRS
gaps{n}(::Type{Bio.Seq.DNAAlphabet{n}})       = DNAGAPS
gaps{n}(::Type{Bio.Seq.RNAAlphabet{n}})       = RNAGAPS
mismatches{n}(::Type{Bio.Seq.DNAAlphabet{n}}) = DNAMISMATCH
mismatches{n}(::Type{Bio.Seq.RNAAlphabet{n}}) = RNAMISMATCH

onehot{I <: Integer}(x::I) = 0x01 << (x-1)

encodeindex{A <: Alphabet}(::Type{A}, x::Bio.Seq.Nucleotide) = (x == DNA_Gap) ? 0 : Bio.Seq.encode(A, x)+1
index{A <: Alphabet}(::Type{A}, func=pairs) = map( x->map(y->encodeindex(A, y), x), func(A) )

Base.reverse( seq::Bio.Seq.ReferenceSequence ) = Bio.Seq.ReferenceSequence( reverse( String( seq ) ) )

type DuplexTrie{A <: Alphabet,K <: Integer}
   fwd::RNATrie{A,K}
   rev::RNATrie{A,K}
   names::Vector{String}
   seqs::Vector{Bio.Seq.Sequence}
   lens::Vector{Int}
   range::UnitRange
   
   function DuplexTrie( seq::Bio.Seq.Sequence, range::UnitRange; seqname="chr" )
      fwd = RNATrie{A,K}( range )
      rev = RNATrie{A,K}( range )
      names = String[seqname]
      seqs  = Bio.Seq.Sequence[seq]
      lens  = Int[length(seq)]
      push!( fwd, seq, one(K) )
      push!( rev, reverse(seq), one(K) )
      return new( fwd, rev, names, seqs, lens, range )
   end
end



function traverse{A}( trie::DuplexTrie{A}, foldrange::UnitRange; 
                      bulge_max::Int=zero(Int), mismatch_max::Int=zero(Int) )

   const pairs_idx      = index(A, pairs)
   const bulges_idx     = index(A, gaps)
   const mismatches_idx = index(A, mismatches)

   const duplex         = RNADuplex()
   const deprange       = trie.range
   const energy_max     = 0.01

   @inline function traverse{A,K}( fwd::TrieNode{A,K}, rev::TrieNode{A,K}, depth::Int64, 
                                   bulge_n::Int64, mismatch_n::Int64,
                                   from_bulge::Bool, bulge_left::Bool )
   
      # recurse through valid watson-crick pairs
      for (l,r) in pairs_idx
         if isa( fwd.next[l], TrieNode{A,K} ) && isa( rev.next[r], TrieNode{A,K} )
            push!( duplex, convert(RNAPair, onehot(l), onehot(r)) ) 
            if depth in deprange && energy(duplex) < trie.range.start*-1
               for (ix,i) in enumerate(fwd.offsets[l]), (jx,j) in enumerate(rev.offsets[r])
                  k = revoffset( j, trie.lens[ rev.metadata[r][jx] ] )
                  if (k - i) + 1 in foldrange
                     #println("$depth: $l + $r @ $i & $k @ $bulge_n @ $mismatch_n && energy=$( energy(duplex))")
                     #println(duplex)
                     
                  end
               end
            end
            traverse( fwd.next[l], rev.next[r], depth + 1,
                       bulge_n, mismatch_n,
                       false, false )
         end
      end
      # recurse through bulges
      if bulge_n < bulge_max && depth > 1 && duplex.energy[end] < energy_max
         for (l,r) in bulges_idx
            if l == 0
               (from_bulge && !bulge_left) && continue # only bulge one way
               if isa( rev.next[r], TrieNode{A,K} )
                  push!( duplex, convert(RNABulge, zero(UInt8), onehot(r)) )
                  traverse( fwd, rev.next[r], depth,
                             bulge_n + 1, mismatch_n,
                             true, true )
               end
            elseif r == 0
               (from_bulge && bulge_left) && continue 
               if isa( fwd.next[l], TrieNode{A,K} )
                  push!( duplex, convert(RNABulge, onehot(l), zero(UInt8)) )
                  traverse( fwd.next[l], rev, depth,
                             bulge_n + 1, mismatch_n,
                             true, false )
               end
            end
         end
      end
      # recurse through mismatches
      if mismatch_n < mismatch_max && depth > 1 && !from_bulge && duplex.energy[end] < energy_max
         for (l,r) in mismatches_idx
            if isa( fwd.next[l], TrieNode{A,K} ) && isa( rev.next[r], TrieNode{A,K} )
             push!( duplex, convert(RNAMismatch, onehot(l), onehot(r)) )
             traverse( fwd.next[l], rev.next[r], depth + 1,
                        bulge_n, mismatch_n + 1,
                        from_bulge, bulge_left )
            end
         end
      end

      pop!( duplex ) # clean up this level
      return
   end

   traverse( trie.fwd.root, trie.rev.root, one(Int),
             zero(Int), zero(Int),
             false, false )
end


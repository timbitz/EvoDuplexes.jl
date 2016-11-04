abstract AbstractTrieNode

typealias Alphabet Bio.Seq.Alphabet

immutable NullTrieNode <: AbstractTrieNode
end

type TrieNode{A<:Alphabet} <: AbstractTrieNode
   next::Vector{AbstractTrieNode}
   offsets::Vector{Vector{Int}}

   function TrieNode()
      next = Vector{AbstractTrieNode}(length(alphabet(A)))
      fill!(next, NullTrieNode())
      return new(next, Vector{Vector{Int}}())
   end
end

type OffsetTrie{A<:Alphabet}
   root::AbstractTrieNode
   range::UnitRange

   OffsetTrie(range::UnitRange) = new(NullTrieNode(), range)
end

function Base.push!{A}( trie::OffsetTrie{A}, seq::Bio.Seq.Sequence )
   if isa( trie.root, NullTrieNode )
      trie.root = TrieNode{A}()
   end
   if length(seq) <= last(trie.range)
      push!( trie.root, seq, trie.range, 1 )
   else
      for i in 1:length(seq)-last(trie.range)+1
         push!( trie.root, seq, trie.range, i )
      end
   end
end

function Base.push!{A}( node::TrieNode{A}, seq::Bio.Seq.Sequence, range::UnitRange, idx::Int, curdepth::Int=1 )
   if curdepth > last(range) || isambiguous( seq[idx] )
      return
   end
   nucidx = Bio.Seq.encode( A, seq[idx] ) + 1
   if isa( node.next[nucidx], NullTrieNode )
      node.next[nucidx] = TrieNode{A}()
   end
   if curdepth in range
      if length(node.offsets) == 0
         node.offsets = Vector{Vector{Int}}(length(alphabet(A)))
      end
      if !isdefined( node.offsets, nucidx )
         node.offsets[nucidx] = Vector{Int}()
      end
      push!( node.offsets[nucidx], idx )
   end
   push!( node.next[nucidx], seq, range, idx + 1, curdepth + 1 )
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

const RNAPAIRS = [(RNA_A, RNA_U), (RNA_U, RNA_A),
                  (RNA_G, RNA_C), (RNA_C, RNA_G),
                  (RNA_G, RNA_U), (RNA_U, RNA_G)]

const RNAGAPS  = [(RNA_A, RNA_Gap), (RNA_Gap, RNA_A),
                  (RNA_C, RNA_Gap), (RNA_Gap, RNA_C),
                  (RNA_G, RNA_Gap), (RNA_Gap, RNA_G),
                  (RNA_U, RNA_Gap), (RNA_Gap, RNA_U)]

typealias PairsType Vector{Tuple{Bio.Seq.Nucleotide, Bio.Seq.Nucleotide}}

pairs{n}(::Type{Bio.Seq.DNAAlphabet{n}}, gap::Bool=false) = gap ? DNAGAPS : DNAPAIRS
pairs{n}(::Type{Bio.Seq.RNAAlphabet{n}}, gap::Bool=false) = gap ? RNAGAPS : RNAPAIRS

encodeindex{A <: Alphabet}(::Type{A}, x::Bio.Seq.Nucleotide) = (x == DNA_Gap) ? 0 : Bio.Seq.encode(A, x)+1
indexpairs{A <: Alphabet}(::Type{A}, gap::Bool=false) = map( x->map(y->encodeindex(A, y), x), pairs(A, gap) )

Base.reverse( seq::Bio.Seq.ReferenceSequence ) = Bio.Seq.ReferenceSequence( reverse( String( seq ) ) )

type DuplexTrie{A <: Alphabet}
   fwd::OffsetTrie{A}
   rev::OffsetTrie{A}
   seq::Bio.Seq.Sequence
   range::UnitRange
   
   function DuplexTrie( seq::Bio.Seq.Sequence, range::UnitRange )
      fwd = OffsetTrie{A}( range )
      rev = OffsetTrie{A}( range )
      push!( fwd, seq )
      push!( rev, reverse(seq) )
      return new( fwd, rev, seq, range )
   end
end

function traverse{A}( trie::DuplexTrie{A}, foldrange::UnitRange, bulge_max=0 )
   traverse( trie.fwd.root, trie.rev.root, 1, trie.range, length(trie.seq), foldrange, 0, bulge_max )
end

function traverse{A}( fwd::TrieNode{A}, rev::TrieNode{A}, 
                      depth::Int, deprange::UnitRange, 
                      len::Int, foldrange::UnitRange, 
                      bulge_n::Int, bulge_max::Int )
   success = false
   # recurse through valid watson-crick pairs
   for (l,r) in indexpairs(A)
      if isa( fwd.next[l], TrieNode{A} ) && isa( rev.next[r], TrieNode{A} )
         success = traverse( fwd.next[l], rev.next[r], 
                             depth + 1, deprange, 
                             len, foldrange,
                             bulge_n, bulge_max )
         if depth in deprange && !success
            for i in fwd.offsets[l], j in rev.offsets[r]
               k = revoffset( j, len )
               if (k - i) + 1 in foldrange
                  println("$depth: $l + $r @ $i & $k @ $bulge_n")
               end
            end
            success = true
         end
      end
   end
   # recurse through bulges
   if bulge_n < bulge_max && depth > 1
      for (l,r) in indexpairs(A, true) # gaps = true
         if l == 0 # gap
            if isa( rev.next[r], TrieNode{A} )
               success = traverse( fwd, rev.next[r], 
                                   depth, deprange, 
                                   len, foldrange,
                                   bulge_n + 1, bulge_max )
            end
         elseif r == 0
            if isa( fwd.next[l], TrieNode{A} )
               success = traverse( fwd.next[l], rev,
                                   depth, deprange,
                                   len, foldrange,
                                   bulge_n + 1, bulge_max )
            end
         end
      end
   end
   # recurse through mismatches
   if mismatch_n < mismatch_max && depth > 1
      #
   end
   success
end


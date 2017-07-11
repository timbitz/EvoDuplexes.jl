
immutable RNAGenomeCoord
   name::String
   offset::Int
   strand::Bool
   length::Int
end

type RNASuffixArray{A<:Bio.Seq.Alphabet,I<:Integer,K<:Integer}
   sai::Vector{I}
   meta::Vector{K}
   depth::Vector{Vector{Bio.Seq.Nucleotide}}
   seqs::Vector{Bio.Seq.Sequence}
   coords::Vector{RNAGenomeCoord}
   length::Int

   function RNASuffixArray( seq::Bio.Seq.Sequence, len::Int, metadata::K=one(K);
                            seqname::String="chr", genomepos::Int=1, strand::Bool=true )
      sai   = SuffixArrays.suffixsort( String(seq) ).index + one(I)
      meta  = K[metadata for i in 1:length(sai)]
      depth = Vector{Vector{Bio.Seq.Nucleotide}}(len)
      for i in 1:len
         depth[i] = map(x->x+i-1 > length(seq) ? gap(A) : seq[x+i-1], sai)
      end
      seqs   = Bio.Seq.Sequence[seq]
      coords = RNAGenomeCoord[RNAGenomeCoord(seqname, genomepos, strand, length(seq))]
      new( sai, meta, depth, seqs, coords, len )
   end
end


function Base.push!{A,I,K}( suf::RNASuffixArray{A,I,K}, seq::Bio.Seq.Sequence )

   # suffix sort
   const len = length(suf.seqs)
   const suftmp  = RNASuffixArray{A,I,K}( seq, suf.length, convert(K, len+1) )

   # initialize
   const newsai    = zeros(I, length(suf.sai) + length(suftmp.sai) )
   const newmeta   = zeros(K, length(suf.sai) + length(suftmp.sai) )
   const newdepth  = Vector{Vector{Bio.Seq.Nucleotide}}(suf.length)
   for i in 1:suf.length
      newdepth[i] = Vector{Bio.Seq.Nucleotide}(length(suf.sai) + length(suftmp.sai))
   end
   const newseqs   = Bio.Seq.Sequence[ suf.seqs; seq ]
   const newcoords = RNAGeomeCoord[ suf.coords; suftmp.coords ]

   # merge
   i,j,k = 1,1,1
   while i <= length(suf.sai) || j <= length(suftmp.sai)
      if j > length(suftmp.sai) || (i <= length(suf.sai) && # does suf have next least suffix?
                                    suf.seqs[ suf.meta[i] ][ suf.sai[i]:end ] < 
                                    suftmp.seqs[ suftmp.meta[j]-len ][ suftmp.sai[j]:end ])
         newsai[k]  = suf.sai[i]
         newmeta[k] = suf.meta[i]
         for d in 1:suf.length
            newdepth[d][k] = suf.depth[d][i]
         end
         i += 1
      else # suftmp has least suffix
         newsai[k]  = suftmp.sai[j]
         newmeta[k] = suftmp.meta[j]
         for d in 1:suf.length
            newdepth[d][k] = suftmp.depth[d][j]
         end
         j += 1
      end
      k += 1
   end

   # set new values of suf
   suf.sai    = newsai
   suf.meta   = newmeta
   suf.depth  = newdepth
   suf.seqs   = newseqs
   suf.coords = newcoords
   suf
end

type RNADuplexArray{A,I,K}
   fwd::RNASuffixArray{A,I,K}
   rev::RNASuffixArray{A,I,K}
   depth::Int

   function RNADuplexArray( seq::Bio.Seq.Sequence, depth::Int;
                            seqname::String="chr", genomepos::Int=1, strand::Bool=true )
      fwd = RNASuffixArray{A,I,K}( seq, depth, seqname=seqname, genomepos=genomepos, strand=strand )
      rev = RNASuffixArray{A,I,K}( reverse(seq), depth, seqname=seqname, genomepos=genomepos, strand=strand )
      return new( fwd, rev, depth )
   end

   function RNADuplexArray( left::Bio.Seq.Sequence, right::Bio.Seq.Sequence, depth::Int;
                            left_seqname::String="chr", right_seqname::String="chr",
                            left_genomepos::Int=1, right_genomepos::Int=1,
                            left_strand::Bool=true, right_strand::Bool=true )
      fwd = RNASuffixArray{A,I,K}( left,  depth, one(K), seqname=left_seqname, genomepos=left_genomepos, strand=left_strand )
      rev = RNASuffixArray{A,I,K}( reverse(right), depth, one(K), seqname=right_seqname, genomepos=right_genomepos, strand=right_strand )
      return new( fwd, rev, depth )
   end

   function RNADuplexArray( seqs::Vector{Bio.Seq.Sequence}, depth::Int;
                            seqnames::Vector{String}=String["chr" for i in 1:length(seqs)],
                            genomepos::Vector{Int}=ones(Int, length(seqs)),
                            strands::Vector{Bool}=trues(length(seqs)))
      # TODO
   end
end

function positions{N <: Bio.Seq.Nucleotide}( nucs::Vector{N}, lo::Int, hi::Int, alpha::Tuple{Vararg{N}} )
   map( x->searchsorted(nucs, x, lo, hi, Base.Order.ForwardOrdering()), alpha ) 
end


function traverse{A,I,K}( dsa::RNADuplexArray{A,I,K}, foldrange::UnitRange;
                          bulge_max::Int=zero(Int), mismatch_max::Int=zero(Int),
                          raw_output::String="", minfold::Float64=-5.0)

   const pairs_idx      = index(A, pairs)
   const bulges_idx     = index(A, gaps)
   const mismatches_idx = index(A, mismatches)

   const alpha          = alphabet(A)

   const duplex         = RNADuplex()
   const deprange       = 1:dsa.depth
   const energy_max     = 0.01

   const intervals      = DuplexCollection{K}()

   if raw_output != ""
      const io = open( raw_output, "w" )
      const stream = ZlibDeflateOutputStream( io )
      const output = true
   else
      const output = false
   end

   @inline function _traverse( fdepth::Int64, rdepth::Int64,
                               frange::UnitRange{Int}, rrange::UnitRange{Int},
                               bulge_n::Int64, mismatch_n::Int64,
                               from_bulge::Bool, bulge_left::Bool,
                               from_mismatch::Bool )

      # obtain position range for current depth
      if fdepth < dsa.fwd.length && rdepth < dsa.rev.length
         const franges = positions( dsa.fwd.depth[fdepth], first(frange), last(frange), alpha )
         const rranges = positions( dsa.rev.depth[rdepth], first(rrange), last(rrange), alpha )
      else
         pop!(duplex)
         return
      end

      # recurse through valid watson-crick pairs
      for (l,r) in pairs_idx
         if length(franges[l]) >= 1 && length(rranges[r]) >= 1
            push!( duplex, convert(RNAPair, onehot(l), onehot(r)) )
            if output
               tab_write( stream, (depth, energy(duplex)))
            end
            if fdepth in deprange && rdepth in deprange && energy(duplex) < minfold
               for (ix,i) in enumerate(dsa.fwd.sai[franges[l]]), (jx,j) in enumerate(dsa.rev.sai[rranges[r]])
                     const fwd_meta = dsa.fwd.meta[ix]
                     const rev_meta = dsa.rev.meta[jx]
                     const len = length( dsa.rev.seqs[ rev_meta ] )
                     const k = revoffset( j, len )
                     fwd_meta == rev_meta && !(((k-rdepth+1) - (i+fdepth-1)) + 1 in foldrange) && continue
                     const newdup = deepcopy(duplex)
                     const fwd_entry = dsa.fwd.coords[ fwd_meta ]
                     const rev_entry = dsa.rev.coords[ rev_meta ]
                     const fg = fwd_entry.offset-1
                     const rg = rev_entry.offset-1
                     push!( intervals, DuplexInterval( Interval(fwd_entry.name, fg+i,   fg+i+fdepth-1, '?', fwd_meta),
                                                       Interval(rev_entry.name, rg+(k-rdepth+1), rg+k, '?', rev_meta),
                                                       newdup ) )
               end
            end
            _traverse(  fdepth + 1, rdepth + 1,
                        franges[l], rranges[r],
                        bulge_n, mismatch_n,
                        false, false, false )
         end
      end

      # recurse through bulges
      if bulge_n < bulge_max && (bulge_n >= 1 ? from_bulge : !from_bulge) &&
         fdepth > 1 && rdepth > 1 && duplex.energy[end] < energy_max
         for (l,r) in bulges_idx
            if l == 0
               (from_bulge && !bulge_left) && continue # only bulge one way
               if length(rranges[r]) >= 1 && rdepth + 1 <= dsa.rev.length
                  push!( duplex, convert(RNABulge, zero(UInt8), onehot(r)) )
                  _traverse( fdepth, rdepth + 1,
                             frange, rranges[r],
                             bulge_n + 1, mismatch_n,
                             true, true, from_mismatch )
               end
            elseif r == 0
               (from_bulge && bulge_left) && continue
               if length(franges[l]) >= 1 && fdepth + 1 <= dsa.fwd.length
                   push!( duplex, convert(RNABulge, onehot(l), zero(UInt8)) )
                   _traverse( fdepth + 1, rdepth,
                              franges[l], rrange,
                              bulge_n + 1, mismatch_n,
                              true, false, from_mismatch )
               end
            end
         end
      end

      # recurse through mismatches
      if mismatch_n < mismatch_max && 
         (mismatch_n >= 1 ? from_mismatch : !from_mismatch) &&
         fdepth > 1 && rdepth > 1 &&
         !from_bulge && duplex.energy[end] < energy_max &&
         fdepth < dsa.fwd.length && rdepth < dsa.rev.length
         for (l,r) in mismatches_idx
            if length(franges[l]) >= 1 && length(franges[r]) >= 1
                push!( duplex, convert(RNAMismatch, onehot(l), onehot(r)) )
                _traverse( fdepth + 1, rdepth + 1,
                           franges[l], rranges[r],
                           bulge_n, mismatch_n + 1,
                           from_bulge, bulge_left, true )
            end
         end
      end

      pop!( duplex ) # clean up this level
      return
   end

   _traverse( one(Int), one(Int),
              1:length(dsa.fwd.sai), 1:length(dsa.rev.sai),
              zero(Int), zero(Int),
              false, false, false )

   if output
      close(stream)
      close(io)
   end

   intervals
end


const SuffixVector = Vector{Vector{Bio.Seq.Nucleotide}}

type RNASuffixArray{A<:Bio.Seq.Alphabet,I<:Integer,K<:Integer}
   sai::Vector{I}
   meta::Vector{K}
   depth::SuffixVector
   species::Vector{SuffixVector}
   length::Int
   isphylo::Bool

   function RNASuffixArray( seq::Bio.Seq.Sequence, len::Int, metadata::K=one(K))
      sai   = SuffixArrays.suffixsort( String(seq) ).index + one(I)
      meta  = K[metadata for i in 1:length(sai)]
      depth = SuffixVector(len)
      for i in 1:len
         depth[i] = map(x->x+i-1 > length(seq) ? gap(A) : seq[x+i-1], sai)
      end
      empty = Vector{SuffixVector}()
      new( sai, meta, depth, empty, len, false )
   end

   function RNASuffixArray( maf::MAFRecord, tree::PhyloTree, len::Int, metadata::K=one(K))
      deletegaps!( maf )
      sai   = SuffixArrays.suffixsort( String(maf[1].sequence) ).index + one(I)
      meta  = K[metadata for i in 1:length(sai)]
      depth = SuffixVector(len)
      for i in 1:len
         depth[i] = map(x->x+i-1 > length(maf[1].sequence) ? gap(A) : maf.species[1].sequence[x+i-1], sai)
      end
      species = Vector{SuffixVector}(length(tree.order)-1)
      m = 2
      for i in 1:length(species)
         curind = tree.index[maf[m].name] - 1
         species[i] = SuffixVector(len)
         if i != curind
            for j in 1:len
               species[i][j] = [gap(A) for i in 1:length(sai)]
            end
         else
           for j in 1:len
              species[i][j] = map(x->x+j-1 > length(maf[m].sequence) ? gap(A) : maf[m].sequence[x+j-1], sai)
           end
           m < length(maf) && (m += 1)
         end
      end
      new( sai, meta, depth, species, len, true )
   end
end

function Base.:(<=)( a::SuffixVector, ai, b::SuffixVector, bi )
   for i in 1:length(a)
      if a[i][ai] < b[i][bi]
         return true
      elseif a[i][ai] > b[i][bi]
         return false
      end
   end
   true
end

function Base.push!{A,I,K}( suf::RNASuffixArray{A,I,K}, seq::Bio.Seq.Sequence; 
                            metadata::K=convert(K,maximum(suf.meta)+1) )

   suf.isphylo && error("Cannot push! single BioSequence when isphylo is enabled!")

   # suffix sort
   const suftmp  = RNASuffixArray{A,I,K}( seq, suf.length, metadata )

   # initialize
   const newsai    = zeros(I, length(suf.sai) + length(suftmp.sai) )
   const newmeta   = zeros(K, length(suf.sai) + length(suftmp.sai) )
   const newdepth  = SuffixVector(suf.length)
   for i in 1:suf.length
      newdepth[i] = Vector{Bio.Seq.Nucleotide}(length(suf.sai) + length(suftmp.sai))
   end

   # merge
   i,j,k = 1,1,1
   while i <= length(suf.sai) || j <= length(suftmp.sai)
      if j > length(suftmp.sai) || (i <= length(suf.sai) && # does suf have next least suffix?
                                    <=( suf.depth, i, suftmp.depth, j ))
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
   suf.sai     = newsai
   suf.meta    = newmeta
   suf.depth   = newdepth
   suf
end


function Base.push!{A,I,K}( suf::RNASuffixArray{A,I,K}, maf::MAFRecord, tree::PhyloTree;
                            metadata::K=convert(K,maximum(suf.meta)+1) )

   suf.isphylo || error("Cannot push! MAF block when isphylo isn't enabled!")
   length(suf.species) == length(tree.order)-1 || error("# species in tree not equal to suffix array!")

   # suffix sort
   const suftmp  = RNASuffixArray{A,I,K}( maf, tree, suf.length, metadata )

   # initialize
   const newsai     = zeros(I, length(suf.sai) + length(suftmp.sai) )
   const newmeta    = zeros(K, length(suf.sai) + length(suftmp.sai) )
   const newdepth   = SuffixVector(suf.length)
   const newspecies = Vector{SuffixVector}(length(tree.order)-1)

   for i in 1:suf.length
      newdepth[i] = Vector{Bio.Seq.Nucleotide}(length(suf.sai) + length(suftmp.sai))
   end
   for i in 1:length(newspecies)
      newspecies[i] = SuffixVector(suf.length)
      for j in 1:suf.length
         newspecies[i][j] = Vector{Bio.Seq.Nucleotide}(length(suf.sai) + length(suftmp.sai))
      end
   end

   # merge
   i,j,k = 1,1,1
   while i <= length(suf.sai) || j <= length(suftmp.sai)
      if j > length(suftmp.sai) || (i <= length(suf.sai) && # does suf have next least suffix?
                                    <=( suf.depth, i, suftmp.depth, j ))
         newsai[k]  = suf.sai[i]
         newmeta[k] = suf.meta[i]
         for d in 1:suf.length
            newdepth[d][k] = suf.depth[d][i]
         end
         for s in 1:length(suf.species)
            for d in 1:suf.length
               newspecies[s][d][k] = suf.species[s][d][i]
            end
         end
         i += 1
      else # suftmp has least suffix
         newsai[k]  = suftmp.sai[j]
         newmeta[k] = suftmp.meta[j]
         for d in 1:suf.length
            newdepth[d][k] = suftmp.depth[d][j]
         end
         for s in 1:length(suf.species)
            for d in 1:suf.length
               newspecies[s][d][k] = suftmp.species[s][d][j]
            end
         end
         j += 1
      end
      k += 1
   end

   # set new values of suf
   suf.sai     = newsai
   suf.meta    = newmeta
   suf.depth   = newdepth
   suf.species = newspecies
   suf
end

immutable RNAGenomeCoord
   name::String
   offset::Int
   strand::Bool
   length::Int
end

type RNADuplexArray{A,I,K}
   fwd::RNASuffixArray{A,I,K}
   rev::RNASuffixArray{A,I,K}
   seqs::Vector{Bio.Seq.Sequence}
   coords::Vector{RNAGenomeCoord}
   depth::Int
   isphylo::Bool

   function RNADuplexArray( seq::Bio.Seq.Sequence, depth::Int;
                            seqname::String="chr", genomepos::Int=1, strand::Bool=true )
      fwd    = RNASuffixArray{A,I,K}( seq, depth )
      rev    = RNASuffixArray{A,I,K}( reverse(seq), depth )
      seqs   = Bio.Seq.Sequence[seq]
      coords = RNAGenomeCoord[RNAGenomeCoord(seqname, genomepos, strand, length(seq))]
      return new( fwd, rev, seqs, coords, depth, false )
   end

   function RNADuplexArray( left::Bio.Seq.Sequence, right::Bio.Seq.Sequence, depth::Int;
                            left_seqname::String="chr", right_seqname::String="chr",
                            left_genomepos::Int=1, right_genomepos::Int=1,
                            left_strand::Bool=true, right_strand::Bool=true )
      fwd    = RNASuffixArray{A,I,K}( left,  depth, convert(K, 1) )
      rev    = RNASuffixArray{A,I,K}( reverse(right), depth, convert(K, 2) )
      seqs   = Bio.Seq.Sequence[left, right]

      coords = RNAGenomeCoord[RNAGenomeCoord(left_seqname, left_genomepos, left_strand, length(left)),
                              RNAGenomeCoord(right_seqname, right_genomepos, right_strand, length(right))]
      return new( fwd, rev, seqs, coords, depth, false )
   end

   function RNADuplexArray( seqs::Vector{Bio.Seq.Sequence}, depth::Int;
                            seqnames::Vector{String}=String["chr" for i in 1:length(seqs)],
                            genomepos::Vector{Int}=ones(Int, length(seqs)),
                            strands::Vector{Bool}=trues(length(seqs)))
      length(seqs) >= 1 || error("Cannot build RNADuplexArray from seq vector of length 0!")
      length(seqs) == length(seqnames) == length(genomepos) == length(strands) || error("all RNADuplexArray input vectors must be equal size!")
      fwd    = RNASuffixArray{A,I,K}( seqs[1], depth )
      rev    = RNASuffixArray{A,I,K}( reverse(seqs[1]), depth )
      coords = RNAGenomeCoord[RNAGenomeCoord(seqnames[i], genomepos[i], strands[i], length(seqs[i])) for i in 1:length(seqs)]
      for i in 2:length(seqs)
         push!( fwd, seqs[i] )
         push!( rev, reverse(seqs[i]) ) 
      end
      return new( fwd, rev, seqs, coords, depth, false )
   end

   function RNADuplexArray( maf::MAFRecord, tree::PhyloTree, depth::Int )
      refmaf = maf[1]
      fwd    = RNASuffixArray{A,I,K}( maf, depth )
      reverse!( maf )
      rev    = RNASuffixArray{A,I,K}( maf, depth )
      seqs   = Bio.Seq.Sequence[refmaf.sequence]
      coords = RNAGenomeCoord[RNAGenomeCoord(refmaf.chr, refmaf.position, refmaf.strand, length(refmaf.sequence))]
      return new( fwd, rev, seqs, coords, depth, true )
   end
end

function Base.push!{A,I,K}( rda::RNADuplexArray{A,I,K}, seq::Bio.Seq.Sequence;
                            seqname::String="chr", genomepos::Int=1, strand::Bool=true )
   rda.isphylo && error("Cannot push single sequence to RNADuplexArray with isphylo=true!")
   push!( rda.fwd, seq )
   push!( rda.rev, reverse(seq) )
   push!( rda.seqs, seq )
   push!( rda.coords, RNAGenomeCoord(seqname, genomepos, strand, length(seq)) )
end

function Base.push!{A,I,K}( rda::RNADuplexArray{A,I,K}, maf::MAFRecord, tree::PhyloTree )
   refmaf = maf[1]
   push!( rda.fwd, maf, tree )
   reverse!( maf )
   push!( rda.rev, maf, tree )
   push!( rda.seqs, refmaf.sequence )
   push!( rda.coords, RNAGenomeCoord(refmaf.chr, refmaf.position, refmaf.strand, length(refmaf.sequence)) )
end

tochar( strand::Bool ) = strand ? '+' : '-'

function positions{N <: Bio.Seq.Nucleotide}( nucs::Vector{N}, lo::Int, hi::Int, alpha::Tuple{Vararg{N}} )
   map( x->searchsorted(nucs, x, lo, hi, Base.Order.ForwardOrdering()), alpha ) 
end


function traverse{A,I,K}( dsa::RNADuplexArray{A,I,K}, foldrange::UnitRange=1:typemax(Int);
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
                     const fwd_entry = dsa.coords[ dsa.fwd.meta[ix] ]
                     const rev_entry = dsa.coords[ dsa.rev.meta[jx] ]
                     const ik = revoffset(i, fwd_entry.length)
                     const jk = revoffset(j, rev_entry.length)
                     const ifirst,ilast = fwd_entry.strand ? (i,i+fdepth-1) : (ik-fdepth+1,ik)
                     const jfirst,jlast = rev_entry.strand ? (jk-rdepth+1,jk) : (j,j+rdepth-1)
                     const fg = fwd_entry.offset-1
                     const rg = rev_entry.offset-1
                     fwd_entry.name == rev_entry.name && !((rg+jfirst - fg+ilast) + 1 in foldrange) && continue
                     if dsa.isphylo
                        const newdup = EvoDuplex(deepcopy(duplex), tree, dsa, ix, jx)
                     else
                        const newdup = deepcopy(duplex)
                     end
                     push!( intervals, DuplexInterval( Interval(fwd_entry.name, fg+ifirst, fg+ilast, tochar(fwd_entry.strand), dsa.fwd.meta[ix]),
                                                       Interval(rev_entry.name, rg+jfirst, rg+jlast, tochar(rev_entry.strand), dsa.rev.meta[jx]),
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

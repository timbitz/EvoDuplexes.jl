
type EvoDuplex <: AbstractDuplex
   duplex::RNADuplex
   alignment::Array{Bio.Seq.Nucleotide,2}
   first::Int
   bracket::Vector{Char}
   score::Float64

   function EvoDuplex( rdup::RNADuplex, tree::PhyloTree, dsa::RNADuplexArray, fwd_index, rev_index )
      alignment = Array{Bio.Seq.Nucleotide,2}(length(tree.order), sum(lengths(rdup.path)))
      first   = 1
      last    = size(alignment, 2)
      bracket = Vector{Char}(size(alignment, 2))
      i,j     = 0,0
      for k in 1:length(rdup.path)
         if isbulge(rdup[k]) && isfiveprime(rdup[k])
            i += 1
            col = Bio.Seq.Nucleotide[dsa.fwd.species[s][i][fwd_index] for s in 1:length(tree.order)-1]
            unshift!(col, dsa.fwd.depth[i][fwd_index])
            alignment[:,first] = col
            bracket[first] = '.'
            first += 1
         elseif isbulge(rdup[k]) && !isfiveprime(rdup[k])
            j += 1
            col = Bio.Seq.Nucleotide[dsa.rev.species[s][j][rev_index] for s in 1:length(tree.order)-1] 
            unshift!(col, dsa.rev.depth[j][rev_index])
            alignment[:,last] = col
            bracket[last] = '.'
            last -= 1
          else
            i += 1
            j += 1
            const lcol = Bio.Seq.Nucleotide[dsa.fwd.species[s][i][fwd_index] for s in 1:length(tree.order)-1]
            const rcol = Bio.Seq.Nucleotide[dsa.rev.species[s][j][rev_index] for s in 1:length(tree.order)-1]
            unshift!(lcol, dsa.fwd.depth[i][fwd_index])
            unshift!(rcol, dsa.rev.depth[j][rev_index])
            #println("dsa.rev.depth $j $rev_index")
            alignment[:,first] = lcol
            alignment[:,last]  = rcol
            bracket[first] = '('
            bracket[last]  = ')'
            first += 1
            last  -= 1
         end
      end
      return new( deepcopy(rdup), alignment, first-1, bracket, 0.0 )
   end

   function EvoDuplex( rdup::RNADuplex, single::PhyloTree, paired::PhyloTree, dsa::RNADuplexArray, fwd_index, rev_index )
      evo = EvoDuplex( rdup, single, dsa, fwd_index, rev_index )
      score!( evo, single, paired )
      return evo
   end
end

path( evo::EvoDuplex )        = evo.duplex.path

energy( evo::EvoDuplex )      = energy( evo.duplex )
energies( evo::EvoDuplex )    = map( energy, [evo.duplex; evo.species] )

npairs( evo::EvoDuplex )      = npairs( evo.duplex )
nmismatches( evo::EvoDuplex ) = nmismatches( evo.duplex )
nbulges( evo::EvoDuplex )     = nbulges( evo.duplex )

strings( evo::EvoDuplex )     = strings( evo.duplex )

Base.push!{NP <: NucleotidePair}( evo::EvoDuplex, pair::NP )         = push!( evo.duplex, pair )
Base.push!{NP <: NucleotidePair}( evo::EvoDuplex, path::Vector{NP} ) = push!( evo.duplex, path )

function join_duplex!( left::EvoDuplex, right::EvoDuplex, npairs, npairs_first, npairs_last )
   push!( left.duplex, path(right.duplex)[npairs+1:end] )
   last  = (size(left.alignment, 2) - left.first) - npairs_last - 1
   first = left.first - npairs_first
   left.alignment = hcat(left.alignment[:,1:first], right.alignment, left.alignment[:,(end-last:end)])
   left.bracket = vcat(left.bracket[1:first], right.bracket, left.bracket[end-last:end])
   left.first += first
   left
end

score( evo::EvoDuplex ) = evo.score

function score!( evo::EvoDuplex, single::PhyloTree, paired::PhyloTree )
   uns_like = 0.0
   str_like = 0.0
   first   = 1
   last    = size(evo.alignment, 2)
   for k in 1:length(evo.duplex.path)
      if isbulge(evo.duplex[k]) && isfiveprime(evo.duplex[k]) 
         like = likelihood( single, evo.alignment[:,first] )
         str_like += log2(like)
         uns_like += log2(like)
         first += 1
      elseif isbulge(evo.duplex[k]) && !isfiveprime(evo.duplex[k])
         like = likelihood( single, evo.alignment[:,last] )
         str_like += log2(like)
         uns_like += log2(like)
         last -= 1
      else
         str_like += log2(likelihood(paired, evo.alignment[:,first], evo.alignment[:,last]))
         uns_like += log2(likelihood(single, evo.alignment[:,first])) + log2(likelihood(single, evo.alignment[:,last]))
         first += 1
         last  -= 1
      end
   end
   evo.score = str_like - uns_like
   evo.score
end

function covariance( evo::EvoDuplex )
   first   = 1
   last    = size(evo.alignment, 2)
   vals    = 0.0
   npairs  = 0
   for k in 1:length(evo.duplex.path)
      if isbulge(evo.duplex[k]) && isfiveprime(evo.duplex[k])
         first += 1
      elseif isbulge(evo.duplex[k]) && !isfiveprime(evo.duplex[k])
         last -= 1
      else
    #     println(evo.alignment[:,first])
    #     println(evo.alignment[:,last])

         cov = covariance_score(evo.alignment[:,first], evo.alignment[:,last])
         #print(" $cov ")
         #maxval = cov > maxval ? cov : maxval
         vals   += cov
         first  += 1
         last   -= 1
         npairs += 1
      end
   end
   #println( " $maxval \n" )
   vals / npairs
end

const PI_MATRIX = [0.0 0.0 0.0 1.0
                   0.0 0.0 1.0 0.0
                   0.0 1.0 0.0 1.0
                   1.0 0.0 1.0 0.0]

delta( a::Bio.Seq.Nucleotide, b::Bio.Seq.Nucleotide ) = a == b ? 1 : 0

function pi_matrix( a::Bio.Seq.Nucleotide, b::Bio.Seq.Nucleotide )
   if isgap( a ) || isgap( b ) 
      return 0.0
   else
      const ind_a = trailing_zeros(reinterpret(UInt8, a)) + 1
      const ind_b = trailing_zeros(reinterpret(UInt8, b)) + 1
      return PI_MATRIX[ ind_a, ind_b ]
   end
end

# Calculated as described by Hofacker et al. Journal of Molecular Biology 2002.
# https://www.ncbi.nlm.nih.gov/pubmed/12079347
function covariance_score( a::Vector{Bio.Seq.Nucleotide}, b::Vector{Bio.Seq.Nucleotide} )
   denom = binomial(length(a), 2)
   numer = 0.0
   pairs = 0
   for i in 1:length(a)-1
      pairs += pi_matrix( a[i], b[i] )
      for j in i+1:length(a)
         numer += (2 - delta(a[i], a[j]) - delta(b[i], b[j])) * pi_matrix( a[i], b[i] ) * pi_matrix( a[j], b[j] )
      end
   end
   pairs += pi_matrix( a[end], b[end] )
   return (numer / denom) * (pairs / length(a))
end

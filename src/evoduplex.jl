
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

function score( evo::EvoDuplex, single::PhyloTree, paired::PhyloTree )
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

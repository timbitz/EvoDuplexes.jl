
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
            println("dsa.rev.depth $j $rev_index")
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

path( evo::EvoDuplex ) = evo.duplex.path

energy( evo::EvoDuplex ) = energy( evo.duplex )
energies( evo::EvoDuplex ) = map( energy, [evo.duplex; evo.species] )

npairs( evo::EvoDuplex )      = npairs( evo.duplex )
nmismatches( evo::EvoDuplex ) = nmismatches( evo.duplex )
nbulges( evo::EvoDuplex )     = nbulges( evo.duplex )

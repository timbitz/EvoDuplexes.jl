
type EvoDuplex <: AbstractDuplex
   species::Vector{RNADuplex}
   index::Vector{UInt16}
   score::Float64

   function EvoDuplex( rdup::RNADuplex, tree::PhyloTree, dsa::RNADuplexArray, fwd_index, rev_index )
      species = RNADuplex[rdup]
      for s in 2:length(tree.order)
         push!(species, RNADuplex())
      end
      i,j = 1,1
      for k in 1:length(rdup)
         const l = reinterpret(UInt8, dsa.fwd.species[1][i])
         const r = reinterpret(UInt8, dsa.rev.species[1][j])
         const lind = l == 0 ? 0 : trailing_zeros(l)+1
         const rind = r == 0 ? 0 : trailing_zeros(r)+1
         if isbulge(rdup.path[i])
            push!( species[2], PAIR_
         else
            
         end
      end
   end
end

energy( pd::EvoDuplex ) = energy( pd.species[1] )
energies( pd::EvoDuplex ) = map( energy, pd.species )


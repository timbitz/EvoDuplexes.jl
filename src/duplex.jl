

# RNADuplex objects represent the path of a duplex as a
# Vector of NucleotidePair's.  Type of the NucleotidePair
# thus defines the operations: Pair, Mismatch, Bulge
# and interior loops are combinations of Mismatches and Bulges
type RNADuplex
   path::Vector{NucleotidePair} # sequence and structure path
   energy::Vector{Float64} # gibbs free energy from nearest-neighbor
   length::Int # shortest length
end

function Base.push!{NP <: NucleotidePair}( duplex::RNADuplex, pair::NP )
   cur_energy = duplex.energy[end]
   if isa( pair, RNAPair )
      if length(duplex.path) >= 1 && isa( duplex.path[end], RNAPair )
         # add stack energy
         cur_energy += stack_energy( pair, duplex.path[end] )
      else
         # add bulge/loop motif from lookup table
      end
   end
   push!(duplex.path, pair)
end



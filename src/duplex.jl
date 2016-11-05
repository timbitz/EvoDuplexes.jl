
# RNADuplex objects represent the path of a duplex as a
# Vector of NucleotidePair's.  Type of the NucleotidePair
# thus defines the operations: Pair, Mismatch, Bulge
# and interior loops are combinations of Mismatches and Bulges
type RNADuplex
   path::Vector{NucleotidePair} # sequence and structure path
   energy::Vector{Float64} # gibbs free energy from nearest-neighbor
   length::Int # shortest length

   RNADuplex() = new(NucleotidePair[], Float64[TURNER_1998_INITIATION], 0)
end

function Base.push!( duplex::RNADuplex, pair::RNAPair )
   cur_energy = duplex.energy[end]
   if length(duplex.path) >= 1 && isa( duplex.path[end], RNAPair )
      # add stack energy
      cur_energy += stack_energy( pair, duplex.path[end] ) # TODO CHECK!
   else
      # add bulge/loop motif from lookup table   
   end
   push!( duplex.energy, cur_energy )
   push!( duplex.path, pair )
end

function Base.push!{NP <: Union{RNAMismatch,RNABulge}}( duplex::RNADuplex, pair::NP )
   push!( duplex.energy, duplex.energy[end] ) # don't calculate energy yet
   push!( duplex.path, pair )
end

function Base.pop!{NP <: NucleotidePair}( duplex::RNADuplex )
   pop!( duplex.path )
   pop!( duplex.energy )
end

# Figure out what kind of bulge/internal loop we have prior to
# the closing RNAPair then score it appropriately
function motif_energy( duplex::RNADuplex, closing_ind::Int )

   bulge_n    = 0
   mismatch_n = 0
   range      = 1:0

   energy   = 0.0 # return value
   # trace back size of internal loop/bulge motif
   for i in (closing_ind-1):-1:1
      if     isa( duplex.path[i], RNAPair )
         range = i:closing_ind
         break
      elseif isa( duplex.path[i], RNAMismatch )
         mismatch_n += 1
      elseif isa( duplex.path[i], RNABulge )
         bulge_n    += 1
      else
         error("$( duplex.path[i] ) is not a valid NucleotidePair!")
      end
   end

   # assign proper energy function to motif
   if     bulge_n > 0 && mismatch_n <= 0
      energy += bulge_energy( duplex, range )
   elseif mismatch_n > 0 && bulge_n <= 0
      energy += internal_energy_symmetric( duplex, range )
   else # mismatch_n > 0 && bulge_n > 0
      energy += internal_energy_asymmetric( duplex, range, mismatch_n, bulge_n )
   end

   energy
end

function bulge_energy( duplex::RNADuplex, range::UnitRange )
   energy = 0.0
   energy += bulge_exception( duplex, range ) # bulge C exception
   if length(range) > 3 # multiple nt bulge
      energy += au_end_penalty( duplex, range )
      energy += bulge_initiation( length(range) - 2 )      
   else # single nt bulge
      energy += bulge_initiation(1)
      energy += stack_energy( duplex.path[range.start], duplex.path[range.stop] )
      energy += # states?
   end
   energy
end

function internal_1x1_energy( duplex::RNADuplex, range::UnitRange )
   (x,y) = split(duplex.path[range.start+1])
   TURNER_2004_INTERNAL_TWO[ index(duplex.path[range.start]),
                             index(duplex.path[range.stop]) ][x,y]
end

function internal_2x2_energy( duplex::RNADuplex, range::UnitRange )
   TURNER_2004_INTERNAL_FOUR[ index(duplex.path[range.start]),
                              index(duplex.path[range.stop]) ][index(range.start+1),index(range.stop-1)]
end

function internal_1x2_energy( duplex::RNADuplex, range::UnitRange )
   # flip?
end

function internal_energy_symmetric( duplex::RNADuplex, range::UnitRange )
   energy = 0.0
   if     length(range) == 3
      return internal_1x1_energy( duplex, range )
   elseif length(range) == 4
      return internal_2x2_energy( duplex, range )
   else
      energy += au_end_penalty( duplex, range, TURNER_2004_INTERNAL_AU_CLOSURE )
      energy += internal_initiation( (length(range) - 2)*2 )
   end
   energy
end

function internal_energy_asymmetric( duplex::RNADuplex, range::UnitRange, mismatch_n, bulge_n )
   energy = 0.0
   if mismatch_n == 1 && bulge_n == 1 # 1x2 lookup table
      energy += internal_1x2_energy( duplex, range )
   else   
      energy += internal_initiation( 2mismatch_n + bulge_n )
      energy += internal_asymmetry_penalty( mismatch_n, bulge_n )
      energy += au_end_penalty( duplex, range, TURNER_2004_INTERNAL_AU_CLOSURE )
      if mismatch_n >= 2
         # add mismatch stabilizers AG/GA GG/UU etc.
      end
   end

   energy
end

function internal_asymmetry_penalty( mismatch_n::Int, bulge_n::Int )
   abs(mismatch_n - (mismatch_n + bulge_n)) * TURNER_2004_ASYMMETRY_PENALTY
end

function bulge_exception( duplex::RNADuplex, range::UnitRange )
   if (duplex.path[range.start] == CG_PAIR && duplex.path[range.start+1] == CB_BULGE) ||
      (duplex.path[range.start] == GC_PAIR && duplex.path[range.start+1] == BC_BULGE) ||
      (duplex.path[range.stop]  == CG_PAIR && duplex.path[range.stop-1]  == CB_BULGE) ||
      (duplex.path[range.stop]  == GC_PAIR && duplex.path[range.stop-1]  == BC_BULGE)
      return TURNER_2004_BULGE_EXCEPTION
   else
      return 0.0
   end
end

function au_end_penalty( duplex::RNADuplex, range::UnitRange, param=TURNER_2004_AU_PENALTY )
   au_end_penalty( duplex.path[range.start], param ) + 
   au_end_penalty( duplex.path[range.stop],  param )
end

function au_end_penalty( x::RNAPair, param=TURNER_2004_AU_PENALTY )
   if x == AU_PAIR || x == UA_PAIR ||
      x == GU_PAIR || x == UG_PAIR
      return param
   else
      return 0.0
   end
end

bulge_initiation( size::Int ) = TURNER_2004_BULGE_INITIATION[ size ]
internal_initiation( size::Int) = TURNER_2004_INTERNAL_INITIATION[ size ]

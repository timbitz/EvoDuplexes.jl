
# RNADuplex objects represent the path of a duplex as a
# Vector of NucleotidePair's.  Type of the NucleotidePair
# thus defines the operations: Pair, Mismatch, Bulge
# and interior loops are combinations of Mismatches and Bulges
type RNADuplex
   path::Vector{NucleotidePair} # sequence and structure path
   energy::Vector{Float64} # gibbs free energy from nearest-neighbor
   length::Int # shortest length

   RNADuplex() = new(NucleotidePair[], Float64[0.0], 0)
end

function Base.push!( duplex::RNADuplex, pair::RNAPair )
   cur_energy = duplex.energy[end]
   if length(duplex.path) == 0
      push!( duplex.path, pair )
      cur_energy += au_end_penalty( duplex.path[1] )
   elseif length(duplex.path) >= 1 && isa( duplex.path[end], RNAPair )
      # add stack energy
      cur_energy += stack_energy( duplex.path[end], pair )
      push!( duplex.path, pair )
   else
      # add bulge/loop motif from lookup table
      push!( duplex.path, pair )
      cur_energy += motif_energy( duplex, length(duplex.path) )
   end
   push!( duplex.energy, cur_energy )
end

function Base.push!{NP <: Union{RNAMismatch,RNABulge}}( duplex::RNADuplex, pair::NP )
   push!( duplex.energy, duplex.energy[end] ) # don't calculate energy yet
   push!( duplex.path, pair )
end

function Base.push!{NP <: NucleotidePair}( duplex::RNADuplex, path::Vector{NP} )
   for i in path
      push!( duplex, i )
   end
end

function Base.pop!( duplex::RNADuplex )
   if length(duplex.path) >= 1
      pop!( duplex.path )
      pop!( duplex.energy )
   end
end

function Base.show( io::IO, duplex::RNADuplex )
   print_index( io, idx ) = print(io, convert(RNANucleotide, UInt8(0x01 << (idx - 1))))
   print(io, "   ")
   for i in duplex.path
      if !isa( i, RNAPair ) && !( isa(i, RNABulge) && !isfiveprime(i) ) 
         char = isa( i, RNAMismatch ) ? split( i )[1] : split( i )
         print_index( io, char )
      else
         print(io, " ")
      end
   end
   print(io, "\n5' ")
   for i in duplex.path
      if isa( i, RNAPair )
         (char,_) = split(i)
         print_index( io, char )
      else
         print(io, " ")
      end
   end
   print(io, " 3'\n3' ")
   for i in duplex.path
      if isa( i, RNAPair )
         (_,char) = split(i)
         print_index( io, char )
      else
         print(io, " ")
      end
   end
   print(io, " 5'  ΔG = ")
   print(io, energy(duplex))
   print(io, "\n   ")
   for i in duplex.path
      if !isa( i, RNAPair ) && !( isa(i, RNABulge) && isfiveprime(i) )
         char = isa( i, RNAMismatch ) ? split( i )[2] : split( i )
         print_index( io, char )
      else
         print(io, " ")
      end
   end
end


@inline function energy( duplex::RNADuplex )
   length(duplex.path) >= 1 ? Float64(signif( duplex.energy[end] + 
                                      au_end_penalty( duplex.path[end] ) +
                                      helix_symmetry( duplex ) + 
                                      TURNER_1998_INITIATION, 5 )) : zero(Float64)
end

# Figure out what kind of bulge/internal loop we have prior to
# the closing RNAPair then score it appropriately
@inline function motif_energy( duplex::RNADuplex, closing_ind::Int )

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

@inline function bulge_energy( duplex::RNADuplex, range::UnitRange )
   energy = 0.0
   energy += bulge_exception( duplex, range ) # bulge C exception
   if length(range) > 3 # multiple nt bulge
      energy += au_end_penalty( duplex, range )
      energy += bulge_initiation( length(range) - 2 )      
   else # single nt bulge
      energy += bulge_initiation(1)
      energy += stack_energy( duplex.path[range.start], duplex.path[range.stop] )
      energy += bulge_states( duplex, range.start+1 )
   end
   energy
end

@inline function internal_1x1_energy( duplex::RNADuplex, range::UnitRange )
   (x,y) = split(duplex.path[range.start+1])
   TURNER_2004_INTERNAL_TWO[ index(duplex.path[range.start]),
                             index(duplex.path[range.stop]) ][x,y]
end

@inline function internal_2x2_energy( duplex::RNADuplex, range::UnitRange )
   TURNER_2004_INTERNAL_FOUR[ index(duplex.path[range.start]),
                              index(duplex.path[range.stop]) ][index(duplex.path[range.start+1]),
                                                               index(duplex.path[range.stop-1])]
end

@inline function internal_1x2_energy( duplex::RNADuplex, range::UnitRange )
   mismatch_pos = isa( duplex.path[range.start+1], RNAMismatch ) ? 1 : 2
   bulge_pos    = isa( duplex.path[range.start+1], RNABulge )    ? 1 : 2
   bulge_five   = isfiveprime( duplex.path[range.start+bulge_pos] )
   @assert mismatch_pos != bulge_pos
   (x_idx, y_idx) = split( duplex.path[range.start+mismatch_pos] )
    a_idx         = split( duplex.path[range.start+bulge_pos] )
   if !bulge_five
      if mismatch_pos < bulge_pos # properly oriented
         return TURNER_2004_INTERNAL_THREE[ index(duplex.path[range.start]),
                                            index(duplex.path[range.stop]), a_idx ][ x_idx, y_idx ]
      else # switch mismatch and bulge
         return TURNER_2004_INTERNAL_THREE[ index(duplex.path[range.start]),
                                            index(duplex.path[range.stop]), y_idx ][ x_idx, a_idx ]
      end
   else
      if mismatch_pos > bulge_pos # properly oriented
         return TURNER_2004_INTERNAL_THREE[ index(flip(duplex.path[range.stop])),
                                            index(flip(duplex.path[range.start])), a_idx ][ y_idx, x_idx ]
      else # switch mismatch and bulge
         return TURNER_2004_INTERNAL_THREE[ index(flip(duplex.path[range.stop])),
                                            index(flip(duplex.path[range.start])), x_idx ][ y_idx, a_idx ]
      end
   end
end

@inline function internal_energy_symmetric( duplex::RNADuplex, range::UnitRange )
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

@inline function internal_energy_asymmetric( duplex::RNADuplex, range::UnitRange, mismatch_n::Int, bulge_n::Int )
   energy = 0.0
   if mismatch_n == 1 && bulge_n == 1 # 1x2 lookup table
      energy += internal_1x2_energy( duplex, range )
   else   
      energy += internal_initiation( 2mismatch_n + bulge_n )
      energy += internal_asymmetry_penalty( mismatch_n, bulge_n )
      energy += au_end_penalty( duplex, range, TURNER_2004_INTERNAL_AU_CLOSURE )
      if mismatch_n >= 2
         energy += internal_energy_mismatch( duplex, range, mismatch_n, bulge_n )
      end
   end

   energy
end

@inline function internal_asymmetry_penalty( mismatch_n::Int, bulge_n::Int )
   abs(mismatch_n - (mismatch_n + bulge_n)) * TURNER_2004_INTERNAL_ASYMMETRY
end

@inline function internal_energy_mismatch( duplex::RNADuplex, range::UnitRange, mismatch_n::Int, bulge_n::Int )
   energy = 0.0
   (i,j)  = five_three( duplex.path, range.start+1:range.stop-1 )
   (h,k)  = five_three( duplex.path, range.start+1:range.stop-1, last=true )
   if mismatch_n == 2 && bulge_n == 1 # 2x3 terminal mismatches
      energy += internal_2x3_stabilizer( duplex.path[range.start], i, j )
      energy += internal_2x3_stabilizer( flip(duplex.path[range.stop]), k, h )
   else 
      energy += internal_other_stabilizer( i, j )
      energy += internal_other_stabilizer( k, h )
   end

   energy
end

@inline function internal_2x3_stabilizer( pair::RNAPair, five::Int, three::Int )
   if is_purine_pyrimidine( pair ) # 5' RY
      if     five == 3 && three == 1 # 5'GA pair
         return -1.2
      end
   else # 5'YR pair
      if     five == 1 && three == 3 # 5'AG pair
         return -0.5
      elseif five == 3 && three == 1 # 5'GA pair
         return -1.1
      end
   end
   if     five == 3 && three == 3 # 5'GG pair
      return -0.8
   elseif five == 4 && three == 4 # 5'UU pair
      return -0.7
   end

   0.0
end

@inline function internal_other_stabilizer( five::Int, three::Int )
   if     five == 1 && three == 3 # 5'AG pair
      return -0.8
   elseif five == 3 && three == 1 # 5'GA pair
      return -1.0
   elseif five == 3 && three == 3 # 5'GG pair
      return -1.2
   elseif five == 4 && three == 4 # 5'UU pair
      return -0.7
   end

   0.0
end

@inline function bulge_exception( duplex::RNADuplex, range::UnitRange )
   if (duplex.path[range.start] == CG_PAIR && duplex.path[range.start+1] == CB_BULGE) ||
      (duplex.path[range.start] == GC_PAIR && duplex.path[range.start+1] == BC_BULGE) ||
      (duplex.path[range.stop]  == CG_PAIR && duplex.path[range.stop-1]  == CB_BULGE) ||
      (duplex.path[range.stop]  == GC_PAIR && duplex.path[range.stop-1]  == BC_BULGE)
      return TURNER_2004_BULGE_EXCEPTION
   else
      return 0.0
   end
end

@inline function bulge_states( duplex::RNADuplex, pos::Int )
   n_states = 1
   isfive = isfiveprime( duplex.path[pos] ) ? 1 : 2
   base   = split(duplex.path[pos])
   i,j    = pos - 1, pos + 1
   while i >= 1 && isa(duplex.path[i], RNAPair ) &&
                 split(duplex.path[i])[isfive] == base
      n_states += 1
      i -= 1
   end
   while j <= length(duplex.path) && isa(duplex.path[j], RNAPair) &&
                                   split(duplex.path[j])[isfive] == base 
      n_states += 1
      j += 1
   end
   RT * log(n_states) * -1
end

@inline function au_end_penalty( duplex::RNADuplex, range::UnitRange, param=TURNER_2004_AU_PENALTY )
   au_end_penalty( duplex.path[range.start], param ) + 
   au_end_penalty( duplex.path[range.stop],  param )
end

@inline function au_end_penalty( x::RNAPair, param=TURNER_2004_AU_PENALTY )
   if x == AU_PAIR || x == UA_PAIR ||
      x == GU_PAIR || x == UG_PAIR
      return convert(Float64, param)
   else
      return 0.0
   end
end

@inline au_end_penalty( x::NucleotidePair ) = 0.0

function helix_symmetry( duplex::RNADuplex )
   i = 1
   j = length(duplex.path)
   while i < j
      if duplex.path[i] != flip(duplex.path[j])
         return 0.0
      end
      i += 1
      j -= 1
   end

   return TURNER_1998_HELICAL_SYMMETRY
end

bulge_initiation( size::Int )   = TURNER_2004_BULGE_INITIATION[ size ]
internal_initiation( size::Int) = TURNER_2004_INTERNAL_INITIATION[ size ]


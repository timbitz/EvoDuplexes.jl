
using Base.Test
using Bio.Seq

include("../src/pairs.jl")
include("../src/energy.jl")
include("../src/duplex.jl")
include("../src/trie.jl")

const pair_set   = [AU_PAIR, UA_PAIR, CG_PAIR, GC_PAIR, GU_PAIR, UG_PAIR]
const mismat_set = [AA_MISMATCH, AG_MISMATCH, AC_MISMATCH, CA_MISMATCH, CC_MISMATCH,
                    CU_MISMATCH, GA_MISMATCH, GG_MISMATCH, UC_MISMATCH, UU_MISMATCH]

@testset "Nucleotide Pairs" begin
   
   @test convert(RNAPair, RNA_A, RNA_U) == AU_PAIR
   @test convert(RNAPair, RNA_U, RNA_A) == UA_PAIR
   @test convert(RNAPair, RNA_C, RNA_G) == CG_PAIR
   @test convert(RNAPair, RNA_G, RNA_C) == GC_PAIR
   @test convert(RNAPair, RNA_G, RNA_U) == GU_PAIR
   @test convert(RNAPair, RNA_U, RNA_G) == UG_PAIR

   @test convert(RNAMismatch, RNA_A, RNA_A) == AA_MISMATCH
   @test convert(RNAMismatch, RNA_A, RNA_G) == AG_MISMATCH
   @test convert(RNAMismatch, RNA_A, RNA_C) == AC_MISMATCH
   @test convert(RNAMismatch, RNA_C, RNA_A) == CA_MISMATCH
   @test convert(RNAMismatch, RNA_C, RNA_C) == CC_MISMATCH
   @test convert(RNAMismatch, RNA_C, RNA_U) == CU_MISMATCH
   @test convert(RNAMismatch, RNA_G, RNA_A) == GA_MISMATCH
   @test convert(RNAMismatch, RNA_G, RNA_G) == GG_MISMATCH
   @test convert(RNAMismatch, RNA_U, RNA_C) == UC_MISMATCH
   @test convert(RNAMismatch, RNA_U, RNA_U) == UU_MISMATCH

   @test convert(RNABulge, RNA_A, RNA_Gap) == AB_BULGE
   @test convert(RNABulge, RNA_C, RNA_Gap) == CB_BULGE
   @test convert(RNABulge, RNA_G, RNA_Gap) == GB_BULGE
   @test convert(RNABulge, RNA_U, RNA_Gap) == UB_BULGE
   @test convert(RNABulge, RNA_Gap, RNA_A) == BA_BULGE
   @test convert(RNABulge, RNA_Gap, RNA_C) == BC_BULGE
   @test convert(RNABulge, RNA_Gap, RNA_G) == BG_BULGE
   @test convert(RNABulge, RNA_Gap, RNA_U) == BU_BULGE

   unique = Dict{Int,Bool}()
   ind_sum    = 0
   for val in [pair_set; mismat_set]
      i = index(val)
      @test !haskey(unique, i)
      unique[i] = true
      ind_sum += i
   end
   @test ind_sum == sum(collect(1:16))

   @test split(AU_PAIR) == (1,4)
   @test split(CG_PAIR) == (2,3)
   @test split(GU_PAIR) == (3,4)

   @test AU_PAIR == flip(UA_PAIR)
   @test UA_PAIR == flip(AU_PAIR)
   @test CG_PAIR == flip(GC_PAIR)
   @test GC_PAIR == flip(CG_PAIR)
   @test UG_PAIR == flip(GU_PAIR)
   @test GU_PAIR == flip(UG_PAIR)

   @test AA_MISMATCH == flip(AA_MISMATCH)
   @test AG_MISMATCH == flip(GA_MISMATCH)
   @test AC_MISMATCH == flip(CA_MISMATCH)
   @test CC_MISMATCH == flip(CC_MISMATCH)
   @test CU_MISMATCH == flip(UC_MISMATCH)
   @test GG_MISMATCH == flip(GG_MISMATCH)
   @test UU_MISMATCH == flip(UU_MISMATCH)

   @test AB_BULGE == flip(BA_BULGE)
   @test CB_BULGE == flip(BC_BULGE)
   @test GB_BULGE == flip(BG_BULGE)
   @test UB_BULGE == flip(BU_BULGE)
   @test BA_BULGE == flip(AB_BULGE)
   @test BC_BULGE == flip(CB_BULGE)
   @test BG_BULGE == flip(GB_BULGE)
   @test BU_BULGE == flip(UB_BULGE)

end

@testset "Energy Parameters" begin
   for i in pair_set, j in pair_set
      if i == GU_PAIR && j == UG_PAIR || i == UG_PAIR && j == GU_PAIR
         @test TURNER_1998_STACK[index(i), index(j)] > 0.0
      else
         @test TURNER_1998_STACK[index(i), index(j)] < 0.0
      end

      @test length(TURNER_2004_INTERNAL_TWO[ index(i), index(j) ]) == 16
      @test sum(TURNER_2004_INTERNAL_TWO[ index(i), index(j) ]) > 0.0
      
      for k in 1:4
         @test length(TURNER_2004_INTERNAL_THREE[ index(i), index(j), k ]) == 16
         @test sum(TURNER_2004_INTERNAL_THREE[ index(i), index(j), k ]) > 0.0
      end

      @test length(TURNER_2004_INTERNAL_FOUR[ index(i), index(j) ]) == 256
      @test sum(TURNER_2004_INTERNAL_FOUR[ index(i), index(j) ]) > 0.0
   end
end

@testset "Duplex Energy Functions" begin

   duplex = RNADuplex()

   # Test simple 1x1 lookup
   duplex.path = [AU_PAIR, AA_MISMATCH, AU_PAIR]
   @test internal_1x1_energy( duplex, 1:3 ) == 1.9
   duplex.path = [UG_PAIR, UU_MISMATCH, UG_PAIR]
   @test internal_1x1_energy( duplex, 1:3 ) == 1.6
   @test motif_energy( duplex, 3 )          == 1.6

   # Test lookup with 1x2 table.
   duplex.path = [AU_PAIR, AA_MISMATCH, BA_BULGE, AU_PAIR]
   @test internal_1x2_energy( duplex, 1:4 ) == 3.7
   # Reverse bulge/mismatch order
   duplex.path = [AU_PAIR, BA_BULGE, AA_MISMATCH, AU_PAIR]
   @test internal_1x2_energy( duplex, 1:4 ) == 3.7
   # Invert bulge side
   duplex.path = [UA_PAIR, AB_BULGE, AA_MISMATCH, UA_PAIR]
   @test internal_1x2_energy( duplex, 1:4 ) == 3.7
   # Invert and reverse
   duplex.path = [UA_PAIR, AA_MISMATCH, AB_BULGE, UA_PAIR]
   @test internal_1x2_energy( duplex, 1:4 ) == 3.7
   # Test range
   duplex.path = [UG_PAIR, UU_MISMATCH, BU_BULGE, UG_PAIR]
   @test internal_1x2_energy( duplex, 1:4 ) == 3.0
   @test motif_energy( duplex, 4 )          == 3.0

   # Test lookup with 2x2 table.
   duplex.path = [AU_PAIR, AA_MISMATCH, AA_MISMATCH, AU_PAIR]
   @test internal_2x2_energy( duplex, 1:4 ) == 2.8
   duplex.path = [UG_PAIR, UU_MISMATCH, UU_MISMATCH, UG_PAIR]
   @test internal_2x2_energy( duplex, 1:4 ) == 3.1
   @test motif_energy( duplex, 4 )          == 3.1

end

@testset "Duplex Calculations" begin
   
   # Duplex examples from http://rna.urmc.rochester.edu/NNDB

   #= Non-self complimentary duplex

     5'GCACG 3'
     3'CGUGC 5'  =#

    duplex = RNADuplex()
    push!( duplex, [GC_PAIR, CG_PAIR, AU_PAIR, CG_PAIR, GC_PAIR] )
    @test length( duplex.path ) == 5
    @test length( duplex.energy ) == length( duplex.path ) + 1
    @test energy( duplex ) == -6.04

   #= Self complimentary duplex

     5'AGCGCU 3'
     3'UCGCGA 5'  =#

    duplex = RNADuplex()
    push!( duplex, [AU_PAIR, GC_PAIR, CG_PAIR, GC_PAIR, CG_PAIR, UA_PAIR] )
    @test length( duplex.path ) == 6
    @test length( duplex.energy ) == length( duplex.path ) + 1
    @test energy( duplex ) == -7.94

   #= 2x2 internal loop

         GA
     5'CA  CG
     3'GU  GC
         AG  =#

    duplex = RNADuplex()
    push!( duplex, [CG_PAIR,AU_PAIR,GA_MISMATCH,AG_MISMATCH,CG_PAIR,GC_PAIR] )
    @test length( duplex.path ) == 6
    @test length( duplex.energy ) == length( duplex.path ) + 1
    @test energy( duplex ) == -1.51
    
   #= 1x5 internal loop

          G
     5'CA     CG
     3'GU     GC
         GAAAG   =#

     duplex = RNADuplex()
     push!( duplex, [CG_PAIR, AU_PAIR, GG_MISMATCH, BA_BULGE, BA_BULGE, BA_BULGE,
                     BG_BULGE, CG_PAIR, GC_PAIR] )
     @test length( duplex.path ) == 9
     @test length( duplex.energy ) == length( duplex.path ) + 1
     @test energy( duplex ) == 4.69
    
   #= 2x3 internal loop with stabilizing mismatches

         GA
     5'CA   CG
     3'GU   GC
         GAG  =#

     duplex = RNADuplex()
     push!( duplex, [CG_PAIR, AU_PAIR, GG_MISMATCH, AA_MISMATCH, BG_BULGE, CG_PAIR, 
                     GC_PAIR] )
     @test length( duplex.path ) == 7
     @test length( duplex.energy ) == length( duplex.path ) + 1
     @test energy( duplex ) == 0.9

   #= Single exception bulge with multiple states

          C
     5'GCC G
     3'CGG C  =#

     duplex = RNADuplex()
     push!( duplex, [GC_PAIR, CG_PAIR, CG_PAIR, CB_BULGE, GC_PAIR] )
     @test length( duplex.path ) == 5
     @test length( duplex.energy ) == length( duplex.path ) + 1
     @test energy( duplex ) == -2.7

   #= Multiple nucleotide bulge

         ACA
     5'GA   G
     3'CU   C  =#

     duplex = RNADuplex()
     push!( duplex, [GC_PAIR, AU_PAIR, AB_BULGE, CB_BULGE, AB_BULGE, GC_PAIR] )
     @test length( duplex.path ) == 6
     @test length( duplex.energy ) == length( duplex.path ) + 1
     @test energy( duplex ) == 5.4

end


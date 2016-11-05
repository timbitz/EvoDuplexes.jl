
using Base.Test

using Bio.Seq
include("../src/pairs.jl")
include("../src/energy.jl")
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


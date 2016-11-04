
using Base.Test
#using RNATrie

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

end

@testset "Energy Parameters" begin
   for i in pair_set, j in pair_set
      if i == GU_PAIR && j == UG_PAIR || i == UG_PAIR && j == GU_PAIR
         @test TURNER_1998_STACK[index(i), index(j)] > 0.0
      else
         @test TURNER_1998_STACK[index(i), index(j)] < 0.0
      end
   end
end



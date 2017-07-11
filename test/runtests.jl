
using Base.Test
using Bio.Seq
using Bio.Intervals
using SuffixArrays

importall Bio.Intervals

include("../src/pairs.jl")
include("../src/energy.jl")
include("../src/duplex.jl")
include("../src/intervals.jl")
include("../src/traverse.jl")
include("../src/trie.jl")
include("../src/suffix.jl")
include("../src/mafreader.jl")

const pair_set   = [AU_PAIR, UA_PAIR, CG_PAIR, GC_PAIR, GU_PAIR, UG_PAIR]
const mismat_set = [AA_MISMATCH, AG_MISMATCH, AC_MISMATCH, CA_MISMATCH, CC_MISMATCH,
                    CU_MISMATCH, GA_MISMATCH, GG_MISMATCH, UC_MISMATCH, UU_MISMATCH]
const bulge_set  = [AB_BULGE, CB_BULGE, GB_BULGE, UB_BULGE, BA_BULGE, BC_BULGE, BG_BULGE, BU_BULGE]

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

   @test is_purine_pyrimidine( AU_PAIR ) == true
   @test is_purine_pyrimidine( GC_PAIR ) == true
   @test is_purine_pyrimidine( GU_PAIR ) == true
   @test is_purine_pyrimidine( UA_PAIR ) == false
   @test is_purine_pyrimidine( CG_PAIR ) == false
   @test is_purine_pyrimidine( UG_PAIR ) == false

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

   @test isfiveprime( AB_BULGE ) == true
   @test isfiveprime( CB_BULGE ) == true
   @test isfiveprime( GB_BULGE ) == true
   @test isfiveprime( UB_BULGE ) == true
   @test isfiveprime( BA_BULGE ) == false
   @test isfiveprime( BC_BULGE ) == false
   @test isfiveprime( BG_BULGE ) == false
   @test isfiveprime( BU_BULGE ) == false

   for b in bulge_set
      @test !ispair(b)
      @test isbulge(b)
      @test !ismismatch(b)
   end
   @test nbulges(bulge_set) == length(bulge_set)
   @test nmismatches(bulge_set) == 0
   @test npairs(bulge_set) == 0

   for m in mismat_set
      @test !ispair(m)
      @test !isbulge(m)
      @test ismismatch(m)
   end
   @test nmismatches(mismat_set) == length(mismat_set)
   @test nbulges(mismat_set) == 0
   @test npairs(mismat_set) == 0

   for p in pair_set
      @test ispair(p)
      @test !isbulge(p)
      @test !ismismatch(p)
   end
   @test npairs(pair_set) == length(pair_set)
   @test nbulges(pair_set) == 0
   @test nmismatches(pair_set) == 0

   @test split(AB_BULGE) == 1
   @test split(CB_BULGE) == 2
   @test split(GB_BULGE) == 3
   @test split(UB_BULGE) == 4
   @test split(BA_BULGE) == 1
   @test split(BC_BULGE) == 2
   @test split(BG_BULGE) == 3
   @test split(BU_BULGE) == 4

   path = [AU_PAIR, CG_PAIR, AB_BULGE]
   @test five_three(path, 1:3) == (1,4)
   @test five_three(path, 1:3, last=true) == (1,3)
   path = [BC_BULGE, CG_PAIR, BA_BULGE]
   @test five_three(path, 1:3) == (2,2)
   @test five_three(path, 1:3, last=true) == (2,1)

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

   # Test helix symmetry
   duplex.path = [AU_PAIR, UA_PAIR, CG_PAIR, GC_PAIR] # no_symmetry
   @test helix_symmetry( duplex ) == 0.0
   duplex.path = [AU_PAIR, CG_PAIR, GC_PAIR, UA_PAIR] # even_symmetry
   @test helix_symmetry( duplex ) == TURNER_1998_HELICAL_SYMMETRY
   duplex.path = [AU_PAIR, CG_PAIR, GU_PAIR, GC_PAIR, UA_PAIR] # odd_symmetry
   @test helix_symmetry( duplex ) == TURNER_1998_HELICAL_SYMMETRY

end

@testset "Duplex Calculations" begin
   
   # Duplex unit-test examples from http://rna.urmc.rochester.edu/NNDB

   #= Non-self complimentary duplex

     5'GCACG 3'
     3'CGUGC 5'  =#

    duplex = RNADuplex()
    push!( duplex, [GC_PAIR, CG_PAIR, AU_PAIR, CG_PAIR, GC_PAIR] )
    @test length( duplex.path ) == 5
    @test length( duplex.energy ) == length( duplex.path ) + 1
    @test energy( duplex ) == -6.01

   #= Self complimentary duplex

     5'AGCGCU 3'
     3'UCGCGA 5'  =#

    duplex = RNADuplex()
    push!( duplex, [AU_PAIR, GC_PAIR, CG_PAIR, GC_PAIR, CG_PAIR, UA_PAIR] )
    @test length( duplex.path ) == 6
    @test length( duplex.energy ) == length( duplex.path ) + 1
    @test energy( duplex ) == -7.98

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
     @test energy( duplex ) == 0.89

   #= Single exception bulge with multiple states

          C
     5'GCC G
     3'CGG C  =#

     duplex = RNADuplex()
     push!( duplex, [GC_PAIR, CG_PAIR, CG_PAIR, CB_BULGE, GC_PAIR] )
     @test length( duplex.path ) == 5
     @test length( duplex.energy ) == length( duplex.path ) + 1
     @test energy( duplex ) == -2.7867

   #= Multiple nucleotide bulge

         ACA
     5'GA   G
     3'CU   C  =#

     duplex = RNADuplex()
     push!( duplex, [GC_PAIR, AU_PAIR, AB_BULGE, CB_BULGE, AB_BULGE, GC_PAIR] )
     @test length( duplex.path ) == 6
     @test length( duplex.energy ) == length( duplex.path ) + 1
     @test energy( duplex ) == 5.34

end

@testset "Duplex Interval Collection" begin

   # Test pushing, overwriting, and ordering
   dup = RNADuplex()
   push!(dup, [CG_PAIR,AU_PAIR,AB_BULGE,GC_PAIR])
   c = DuplexCollection{String}()
   i = DuplexInterval(Interval("c",1,6,'+',"genea"), 
                      Interval("c",50,55,'+',"genea"), dup)
   push!( c, i )
   @test length(c.names["c"]) == 1
   @test length(c.names["c"][one(UInt64)]) == 1

   dup2 = deepcopy( dup )
   push!( dup2, CG_PAIR )
   @test energy(dup2) < energy(dup)
   i2 = DuplexInterval(Interval("c",2,7,'+',"genea"), 
                       Interval("c",20,25,'+',"genea"), dup2)
   push!( c, i2 )
   @test length(c.names["c"]) == 1
   @test length(c.names["c"][one(UInt64)]) == 2

   i3 = DuplexInterval(Interval("c",1,7,'+',"genea"), 
                       Interval("c",49,55,'+',"genea"), dup2)
   push!( c, i3 )
   @test length(c.names["c"]) == 1
   @test length(c.names["c"][one(UInt64)]) == 2

   @test isordered(c)

   # Test addition-of-larger overwriting duplex
   c = DuplexCollection{String}()
   dup = RNADuplex()
   push!(dup, [CG_PAIR,GC_PAIR,CG_PAIR])
   i = DuplexInterval(Interval("c",1,5,'+',"genea"),
                      Interval("c",100,105,'+',"genea"), dup)
   dup2 = deepcopy(dup)
   push!(dup2, [GC_PAIR,CG_PAIR])
   @test energy(dup2) < energy(dup)
   i2 = DuplexInterval(Interval("c",6,10,'+',"genea"),
                       Interval("c",106,110,'+',"genea"), dup2)
   push!( c, i  )
   push!( c, i2 )
   @test length(c.names["c"]) == 1
   @test length(c.names["c"][one(UInt64)]) == 2

   worse = RNADuplex()
   push!(worse, [GU_PAIR,AB_BULGE,GB_BULGE,UA_PAIR,UG_PAIR,CB_BULGE,UA_PAIR])
   @test energy(worse) > energy(dup)
   i3 = DuplexInterval(Interval("c",1,10,'+',"genea"),
                       Interval("c",100,110,'+',"genea"), worse)

   push!( c, i3 ) # less favorable energy should not get pushed.
   @test length(c.names["c"]) == 1
   @test length(c.names["c"][one(UInt64)]) == 2
   firstint = c.names["c"][one(UInt64)][1]
   @test firstint.first.first == 1 && firstint.first.last == 5

   better = RNADuplex()
   push!(better, [CG_PAIR,GC_PAIR,CG_PAIR,CG_PAIR,GC_PAIR,CG_PAIR])
   @test energy(better) < energy(dup2)
   i3 = DuplexInterval(Interval("c",1,10,'+',"genea"),
                       Interval("c",100,110,'+',"genea"), better)
   push!( c, i3 ) # better energy should overwrite overlapping
   @test length(c.names["c"]) == 1
   firstint = c.names["c"][one(UInt64)][1]
   @test firstint.first.first == 1 && firstint.first.last == 10
   @test length(c.names["c"][one(UInt64)]) == 1

   # Test addition of larger overwriting duplex, and subsequent deletion because of smaller higher energy duplex
   c = DuplexCollection{String}()
   push!( c, i  )
   push!( c, i2 )
   middle = RNADuplex()
   push!(middle, [CG_PAIR,GC_PAIR,CG_PAIR, GC_PAIR])
   @test energy(dup2) < energy(middle) < energy(dup)
   i3 = DuplexInterval(Interval("c",1,10,'+',"genea"),
                       Interval("c",100,110,'+',"genea"), middle)
   push!( c, i3 )
   @test length(c.names["c"]) == 1
   firstint = c.names["c"][one(UInt64)][1]
   @test firstint.first.first == 6 && firstint.first.last == 10
   @test length(c.names["c"][one(UInt64)]) == 1

end

@testset "Duplex Intervals Stitching" begin

   # test duplex strings function

   c = DuplexCollection{String}()
   dup = RNADuplex()
   push!(dup, [CG_PAIR,GC_PAIR,CG_PAIR,CG_PAIR,GC_PAIR,GC_PAIR])
   i = DuplexInterval(Interval("c",  1,  6,'+',"genea"),
                      Interval("c",105,110,'+',"genea"), dup)
   dup2 = RNADuplex()
   push!(dup2, [GC_PAIR,GC_PAIR,UG_PAIR,AU_PAIR,CG_PAIR,GC_PAIR])
   i2 = DuplexInterval(Interval("c",  5, 10,'+',"genea"),
                       Interval("c",101,106,'+',"genea"), dup2)

   @test !isnull( stitch(i, i2, 3, 3) )
   @test !isnull( stitch(i2, i, 3, 3) )   

   # test stitching with bulges
   # test negative stitching cases
end

@testset "RNA Trie Building" begin

   # Vestigial tests for a deprecated data structure

   A = Bio.Seq.DNAAlphabet{2}
   trie = RNATrie{A,String}( 1:3 )
   @test trie.range == 1:3
   @test isa( trie.root, NullTrieNode ) == true
   push!( trie, dna"ACGT", "DNA" )
   @test isa( trie.root, TrieNode{A} )  == true
   for i in 1:4
      @test trie.root.offsets[i] == [i]
   end
   @test trie.root.next[1].offsets[2] == [2]
   @test trie.root.next[1].next[2].offsets[3] == [3]
   for i in 1:4
      @test isa( trie.root.next[1].next[2].next[3].next[i], NullTrieNode ) == true
   end
   @test trie.root.next[2].next[3].offsets[4] == [4]
   for i in 1:4
      @test isa( trie.root.next[2].next[3].next[4].next[i], NullTrieNode ) == true
   end
   @test nodecount( trie ) == 9

   push!( trie, ReferenceSequence("ACGT"), "REF" )
   @test isa( trie.root, TrieNode{A} )  == true
   for i in 1:4
      @test trie.root.offsets[i] == [i,i]
   end
   @test trie.root.next[1].offsets[2] == [2,2]
   @test trie.root.next[1].next[2].offsets[3] == [3,3]
   for i in 1:4
      @test isa( trie.root.next[1].next[2].next[3].next[i], NullTrieNode ) == true
   end
   @test trie.root.next[2].next[3].offsets[4] == [4,4]
   for i in 1:4
      @test isa( trie.root.next[2].next[3].next[4].next[i], NullTrieNode ) == true
   end 
   @test trie.root.next[2].next[3].metadata[4] == String["DNA","REF"]
   @test nodecount( trie ) == 9
end

@testset "Duplex Trie Building and Traversal" begin

   # Duplex Trie deprecated...

   seq = dna"AAATGATGCCGCAGGGGGGGGGGTGCGGCAATCATTT"
   trie = DuplexTrie{DNAAlphabet{2},UInt8}( seq, 8:16 )
   @test length(traverse( trie, 8:100, bulge_max=0 )) == 0
   val = traverse( trie, 8:100, bulge_max=1 )
   @test length(val) == 1

end

@testset "RNASuffixArray Building and Traversal" begin

   seq = dna"AAATGATGCCGCAGGGGGGGGGGTGCGGCAATCATTT"
   rsa = DuplexTrie{DNAAlphabet{2},UInt8,UInt8}( seq, 25 )
   @time length(traverse( rsa, 8:50, bulge_max=0 )) == 0
   @time length(traverse( rsa, 8:50, bulge_max=1 )) == 1

   # test fwd and rev sequences for intermolecular constructor
   # test vector of sequences constructor

   # test larger scale example hnRNP D and genome offsets

   # test recursive stitching
   # test iterval collection and filtering
end

@testset "MAF Parser" begin
   mafheader = """
##maf version=1 scoring=tba.v8 
 # tba.v8 (((human chimp) baboon) (mouse rat)) 
 # multiz.v7
 # maf_project.v5 _tba_right.maf3 mouse _tba_C
 # single_cov2.v4 single_cov2 /dev/stdin
"""
   maflines = """
a score=23262.0
s hg16.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

a score=5062.0
s hg16.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

a score=6636.0
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA
"""
   
   maffile = mafheader * "\n" * maflines

   # Single buffer:

   headreader = MAFReader(BufferedInputStream(IOBuffer(maffile)))
   noheadreader = MAFReader(BufferedInputStream(IOBuffer(maflines)))

   discard_header!(headreader)
   @test headreader.stream.buffer[1:10] == noheadreader.stream.buffer[1:10]   

   headreader = MAFReader(BufferedInputStream(IOBuffer(maffile)))

   while !done(headreader)
      headrec   = MAFRecord()
      noheadrec = MAFRecord()
      read!(headreader, headrec)
      read!(noheadreader, noheadrec)
      @test length(noheadrec) == length(headrec)
      @test length(noheadrec) > 1
      for i in 1:length(headrec)
         @test headrec[i] == noheadrec[i]
      end
   end
   @test done(headreader)
   @test done(noheadreader)

   # Multi buffer

   longreader = MAFReader(BufferedInputStream(IOBuffer(maffile)))
   shortreader = MAFReader(BufferedInputStream(IOBuffer(maffile), 512))
   
   while !done(longreader) 
      longrec  = MAFRecord()
      shortrec = MAFRecord()
      read!(longreader, longrec)
      read!(shortreader, shortrec)
      @test length(shortrec) == length(longrec)
      @test length(shortrec) > 1
      for i in 1:length(shortrec)
         @test shortrec[i] == longrec[i]
      end
   end
   @test done(shortreader)
   @test done(longreader)

   # Bad Records
   maffile = """
##maf version=1 scoring=roast.v3.3
a score=781.000000
s hg19.chr22                16092278 1 +  51304566 C
s panTro4.chr22             14470458 1 +  49737984 C
s nomLeu3.chrUn_GL397501_2   1216368 1 -   1485088 C
s chlSab1.chrUn_KE147573       66869 1 -     85383 C
s monDom5.chr4              96780965 1 - 435153693 T
s danRer7.chr3                323838 1 -  63268876 C
s petMar2.GL501154             43910 1 -     71231 t

a score=-50707.000000
s hg19.chr22                     16092279 32 +  51304566 AGTCCCCAGGTGGT------CATGACACCTCAATTGGA
s panTro4.chr22                  14470459 32 +  49737984 AGTCCCCAGGTGGT------CATGACACCTCAATTGGA
s nomLeu3.chrUn_GL397501_2        1216369 32 -   1485088 AGTCCCCAGGTGGT------CATGACACTTCAACTGGA
s chlSab1.chrUn_KE147573            66870 32 -     85383 AGTCCCCAGGTGAT------CATGACACCTTAACTGGA
s monDom5.chr4                   96780966 24 - 435153693 TGACTCCAGCCAGA------GACGGAGCCC--------
s falChe1.KB397375                  11884 29 +     36355 GACCCTTGTGAAAA------CATGAGGGCTGTGAA---
s falPer1.KB390849                  40484 13 +     49439 AGCCCTTGTCTGA-------------------------
s ficAlb2.chrUn_KE165385           131290 31 -    141484 -cctccccgagtat------cCTGCGACCCCAGATCGA
s geoFor1.JH741012                   5677 32 -     10893 GGGCCCTTTGGGAA------GAGGACACTGCGACTGGA
s taeGut2.chr26_EQ833276_random      8281 32 -     85973 AGCCCCTCAGCTGT------CCCCACACCTCTCCTGCA
s pseHum1.KB221494                  45158 24 +     77128 gaccccacaatgac-------------cctcggttga-
s galGal4.chrUn_JH375878              275 32 -      4229 AGCCCAGTggggga------ccccggagccaaagggga
s pelSin1.JH208549                  94224 38 +    588315 GGCCCCAGGGCCAAATCCGTCCCTGGGTCCATATTCAA
s danRer7.chr3                     323839 32 -  63268876 ATTATCCAACTAGA------CTTGGCTGTTTAGAAGAC
s petMar2.GL501154                  43911 30 -     71231 --agttcacgtgag------ctctacaactcaggtgaa
"""

   noshort = MAFReader(BufferedInputStream(IOBuffer(maffile)))
   noshort.minlen = 5
   noshortrec = MAFRecord()
   read!(noshort, noshortrec)
   @test done(noshort)
   @test length(noshortrec) == 15

   noweak = MAFReader(BufferedInputStream(IOBuffer(maffile)))
   noweak.minscore = 0.0
   noweakrec = MAFRecord()
   read!(noweak, noweakrec)
   @test !done(noweak)
   @test length(noweakrec) == 7
   res = read!(noweak, MAFRecord())
   @test done(noweak)
   @test length(res) == 0

   #Trailing newline
   withnew = MAFReader(BufferedInputStream(IOBuffer(maffile * "\n")))
   withnew.minlen = 5
   withnewrec = MAFRecord()
   read!(withnew, withnewrec)   
   @test done(withnew)
   @test length(withnewrec) == 15

end


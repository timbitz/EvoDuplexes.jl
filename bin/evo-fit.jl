#!/usr/bin/env julia

const dir = abspath( splitdir(@__FILE__)[1] )

tic()
println( STDERR, "EvoDuplexes.jl loading and compiling... " )

using ArgParse
using BufferedStreams
using Libz
using Bio.Intervals
using Bio.Seq

unshift!( LOAD_PATH, dir * "/../src" )
using EvoDuplexes

function parse_cmd()
   s = ArgParseSettings()
   @add_arg_table s begin
      "--tree"
       help     = "Phylogenetic tree with neutral branch lengths (in newick format)"
       arg_type = String
       required = true
      "--regions"
       help     = "BED file containing regions for training"
       arg_type = String
       required = true
      "--maf"
       help     = "Directory with MAF files named by chromosome (chr1.maf.gz...)"
       arg_type = String
       required = true
       default  = "../maf"
   end
   return parse_args(s)
end


function main()
   args = parse_cmd()

   println(STDERR, " $( round( toq(), 6 ) ) seconds" )

   neutraltree = parsenewick(chomp(String(read(open(args["tree"])))))
   pairtree = deepcopy(neutraltree)
   constree = deepcopy(neutraltree)

   extend_branches!( constree,    0.5 )
   extend_branches!( neutraltree, 1.0 )
   extend_branches!( pairtree,    0.5 )

   set_prob_mat!( neutraltree, GTR_SINGLE_Q )
   set_prob_mat!( constree,    GTR_SINGLE_Q )
   set_prob_mat!( pairtree,    GTR_PAIRED_Q )

   regs = loadbed( args["regions"] )

   out = BufferedOutputStream(open("training.jlt", "w"))

   for chr in keys(regs.trees)
      println("Reading $(args["maf"])/$chr.maf.gz")
      reader = ZlibInflateInputStream(open("$(args["maf"])/$chr.maf.gz"), bufsize=129000) |> MAFReader
      @time maf = readmaf!( reader, neutraltree.index, minlength=50, minspecies=10, minscore=0.0, regionbool=true, regions=regs )
      for r in collect(regs)
         haskey( maf.trees, r.seqname ) || continue
         col = collect(intersect( maf, r ))
         length(col) > 1 || continue
         if r.strand == STRAND_NEG
            for i in col
               reverse_complement!( i.metadata )
            end
         end
         println("Building DuplexArray for $(r.metadata.name)")
         rda = RNADuplexArray{DNAAlphabet{2},UInt16,UInt16}( col[1].metadata, neutraltree, 50 )
         for i in 2:length(col)
            push!( rda, col[i].metadata, neutraltree )
         end

         @time trav = stitch(traverse( rda, 3:10000, single=neutraltree, bulge_max=3, mismatch_max=3, minfold=-4.0 ))
         for i in collect(trav)
            const conscore = (score( i.duplex, constree, gapdenom=1.1 ) - score( i.duplex, neutraltree )) / npairs( i.duplex )
            score!( i.duplex, neutraltree, pairtree )
            const neutral = score( i.duplex ) / npairs( i.duplex )
            score!( i.duplex, constree, pairtree )
            const slow    = score( i.duplex ) / npairs( i.duplex )
            const ent     = dinucleotide_entropy(i.duplex.duplex)
            const meancov = covariance( i.duplex ) * 10
            const len     = npairs( i.duplex )
            const dist    = distance( i )
            const deltag  = energy( i.duplex )
            write( out, string(conscore) ); write( out, '\t' )
            write( out, string(neutral) ); write( out, '\t' )
            write( out, string(slow) ); write( out, '\t' )
            write( out, string(deltag) ); write( out, '\t' )
            write( out, string(len) ); write( out, '\t' )
            write( out, string(ent) ); write( out, '\t' )
            write( out, string(meancov) ); write( out, '\t' )
            write( out, string(dist) ); write( out, '\n' )
         end
      end
   end
      
   close(out)
end

main()


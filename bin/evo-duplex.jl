#!/usr/bin/env julia

const dir = abspath( splitdir(@__FILE__)[1] )

tic()
println( STDERR, "EvoDuplexes.jl loading and compiling... " )

using ArgParse
using BufferedStreams
using Libz
using GenomicFeatures
using BioSequences
using DataFrames
using Gadfly
using Distances
using Clustering
using MultivariateStats
using StatsBase

unshift!( LOAD_PATH, dir * "/../src" )
using EvoDuplexes
using EvoDuplexes: score, score!

function parse_cmd()
   s = ArgParseSettings()
   @add_arg_table s begin
      "--tree"
       help     = "Phylogenetic tree with neutral branch lengths (in newick format)"
       arg_type = String
       required = true
      "--cons-regions"
       help     = "BED file containing conserved regions"
       arg_type = String
       required = false
      "--gene-regions"
       help     = "BED file containing gene units to allow long-range folding within"
       arg_type = String
       required = true
      "--maf"
       help     = "Directory with MAF files named by chromosome (chr1.maf.gz...)"
       arg_type = String
       required = true
       default  = "../maf"
      "--model-load"
       help     = "Load pre-trained IsolationForest models, .evt.jls"
       arg_type = String
       required = false
      "--model-data"
       help     = "Load training data from `.jlt` file, output .evt.jls file"
       arg_type = String
       required = false
      "--model-train"
       help     = "Train the IsolationForest models from scratch using input regions, output .evt.jls file"
       action   = :store_true
      "--output","-o"
       help     = "Output prefix"
       arg_type = String
       default  = "evoduplex"
      "--cons-branch"
       help     = "When conserved regions are given, use this branch length multiplier"
       arg_type = Float64
       default  = 0.33
    end
   return parse_args(s)
end


function main()

   istrue( arg ) = (arg != nothing)

   args = parse_cmd()

   println(STDERR, " $( round( toq(), 6 ) ) seconds" )

   neutraltree = parsenewick(String(chomp(String(read(open(args["tree"]))))))
   pairtree = deepcopy(neutraltree)
   constree = deepcopy(neutraltree)

   extend_branches!( neutraltree, 1.0 )
   extend_branches!( constree,    args["cons-branch"] )
   extend_branches!( pairtree,    args["cons-branch"] )

   set_prob_mat!( neutraltree, GTR_SINGLE_Q )
   set_prob_mat!( constree,    GTR_SINGLE_Q )
   set_prob_mat!( pairtree,    GTR_PAIRED_Q )

   secbool   = false
   conserved = false
   if istrue( args["cons-regions"] )
      regs  = loadbed( args["cons-regions"], expandfirst=25, expandlast=25 )
      conserved = true
      if istrue( args["gene-regions"] )
         genes   = loadbed( args["gene-regions"], expandfirst=2500, expandlast=2500 )
         secbool = true
      else
         error("No specific `--gene-regions` given! Can't fold entire chromosomes..")
      end
   elseif istrue( args["gene-regions"] )
      genes = loadbed( args["gene-regions"] )
      regs  = genes
   else
      error("No specific regions given! Can't fold entire chromosomes..")
   end

   train = false
   if istrue( args["model-load"] )
      println("Loading model from $(args["model-load"]) file...")
      model = open(deserialize, args["model-load"])
   elseif istrue( args["model-data"] )
      println("Training model from .jlt datafile $(args["model-data"])...")
      tabdata = open(readdlm, args["model-data"])
      model = DistanceForest()
      train!( model, tabdata )
   elseif args["model-train"]
      train = true
   else
      error("No --model-file or --model-data given!, if training is desired use --model-train!")
   end

   datout = BufferedOutputStream(open("$(args["output"])-raw.jlt", "w"))
   bedout = open("$(args["output"]).bed", "w")
   duplexes = Vector{DuplexInterval}() 

   for chr in keys(genes.trees)
      println("Reading $(args["maf"])/$chr.maf.gz")
      reader = ZlibInflateInputStream(open("$(args["maf"])/$chr.maf.gz"), bufsize=129000) |> MAFReader
      @time maf = readmaf!( reader, neutraltree.index, minlength=30, minspecies=2, minscore=-Inf, regionbool=true, regions=genes, secbool=secbool, secregions=regs )

      for r in collect(genes)
         haskey( maf.trees, r.seqname ) || continue
         col = collect(intersect( maf.trees[r.seqname], r ))
         length(col) >= 1 || continue
         if r.strand == STRAND_NEG
            for i in col
               EvoDuplexes.reverse_complement!( i.metadata )
            end
         end
         println("Building DuplexArray for $(r.metadata)")
         rda = RNADuplexArray{DNAAlphabet{2},UInt16,UInt16}( col[1].metadata, neutraltree, 50 )
         for i in 2:length(col)
            push!( rda, col[i].metadata, neutraltree )
         end

         @time trav = collect(traverse( rda, 1:5000, single=neutraltree, bulge_max=3, mismatch_max=3, minfold=-8.0 ))

         println("Calculating structural conservation across the phylogeny...")
         for i in trav
            #const conscore = (score( i.duplex, constree, gapdenom=1.25 ) - score( i.duplex, neutraltree )) / npairs( i.duplex )
            if !conserved
               score!( i.duplex, neutraltree, pairtree )
               const cons = score( i.duplex ) / npairs( i.duplex )
            else
               score!( i.duplex, constree, pairtree )
               const cons = score( i.duplex ) / npairs( i.duplex )
            end
            const ent     = dinucleotide_entropy(i.duplex.duplex) / (npairs( i.duplex ) + 5)
            const meancov = covariance( i.duplex ) * 10
            const len     = npairs( i.duplex )
            const dist    = distance( i )
            const deltag  = energy( i.duplex )

            write( datout, string(dist) ); write( datout, '\t' )
            write( datout, string(cons) ); write( datout, '\t' )
            write( datout, string(deltag) ); write( datout, '\t' )
            write( datout, string(len) ); write( datout, '\t' )
            write( datout, string(ent) ); write( datout, '\t' )
            write( datout, string(meancov) ); write( datout, '\n' ) 
         end

         duplexes = vcat(duplexes, trav)

      end
   end
    
   close(datout)
   tabdata = open(readdlm, "$(args["output"])-raw.jlt")

   if train
      println("Training IsolationForest models...")
      model = DistanceForest()
      train!( model, tabdata )
   end

   println("Identifying statistical anomalies...")
   sig = EvoDuplexes.predict( model, tabdata )

   println("Found $(length(sig)) significant duplexes..")

   for i in sig
      writebed( bedout, duplexes[i], "EvoDuplex" )
   end

   close(bedout)
end


@time main()


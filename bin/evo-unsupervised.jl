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
      "--regions"
       help     = "BED file containing regions for training"
       arg_type = String
       required = true
      "--maf"
       help     = "Directory with MAF files named by chromosome (chr1.maf.gz...)"
       arg_type = String
       required = true
       default  = "../maf"
      "--output","-o"
       help     = "Output prefix"
       arg_type = String
       default  = "evoduplex"
      "--train"
       help     = "Train the Distance-specific IsolationForest models"
       action   = :store_true
   end
   return parse_args(s)
end


function main()
   args = parse_cmd()

   println(STDERR, " $( round( toq(), 6 ) ) seconds" )

   neutraltree = parsenewick(String(chomp(String(read(open(args["tree"]))))))
   pairtree = deepcopy(neutraltree)
   constree = deepcopy(neutraltree)

   extend_branches!( constree,    0.33 )
   extend_branches!( neutraltree, 1.0 )
   extend_branches!( pairtree,    0.33 )

   set_prob_mat!( neutraltree, GTR_SINGLE_Q )
   set_prob_mat!( constree,    GTR_SINGLE_Q )
   set_prob_mat!( pairtree,    GTR_PAIRED_Q )

   regs = loadbed( args["regions"], expandfirst=25, expandlast=25 )

   datout = BufferedOutputStream(open("$(args["output"])-raw.jlt", "w"))
   bedout = open("$(args["output"]).bed", "w")

   for chr in keys(regs.trees)
      println("Reading $(args["maf"])/$chr.maf.gz")
      reader = ZlibInflateInputStream(open("$(args["maf"])/$chr.maf.gz"), bufsize=129000) |> MAFReader
      @time maf = readmaf!( reader, neutraltree.index, minlength=30, minspecies=2, minscore=-Inf, regionbool=true, regions=regs )
      println(maf)
      for r in collect(regs)
         haskey( maf.trees, r.seqname ) || continue
         println("R:"); println(r)
         col = collect(intersect( maf.trees[r.seqname], r ))
         println("COL:"); println(col)
         length(col) >= 1 || continue
         if r.strand == STRAND_NEG
            for i in col
               EvoDuplexes.reverse_complement!( i.metadata )
            end
         end
         println("IinCOL:")
         for i in col
            println(i.metadata.species)
         end
         println("Building DuplexArray for $(r.metadata)")
         rda = RNADuplexArray{DNAAlphabet{2},UInt16,UInt16}( col[1].metadata, neutraltree, 50 )
         for i in 2:length(col)
            push!( rda, col[i].metadata, neutraltree )
         end

         @time trav = collect(traverse( rda, 1:5000, single=neutraltree, bulge_max=3, mismatch_max=3, minfold=-10.0 ))
         sort!( trav, lt=(a,b)->(energy(a.duplex) < energy(b.duplex)) )
         distnum = min(10000, length(trav))
         halfnum = distnum >> 1
         indices = collect(1:halfnum)
         indices = [indices; sample(collect(halfnum+1:distnum), distnum - halfnum, replace=false)]

         dcons = zeros(length(indices))
         dslow = zeros(length(indices))
         ddelg = zeros(length(indices))
         ddist = zeros(length(indices))
         dentr = zeros(length(indices))

         for k in 1:length(indices)
            j = indices[k]
            i = trav[j]
            const conscore = (score( i.duplex, constree, gapdenom=1.25 ) - score( i.duplex, neutraltree )) / npairs( i.duplex )
            score!( i.duplex, neutraltree, pairtree )
            const neutral = score( i.duplex ) / npairs( i.duplex )
            score!( i.duplex, constree, pairtree )
            const slow    = score( i.duplex ) / npairs( i.duplex )
            const ent     = dinucleotide_entropy(i.duplex.duplex) / (npairs( i.duplex ) + 5)
            const meancov = covariance( i.duplex ) * 10
            const len     = npairs( i.duplex )
            const dist    = distance( i )
            const deltag  = energy( i.duplex )

            dcons[k] = conscore 
            dslow[k] = slow
            ddelg[k] = deltag 
            ddist[k] = dist 
            dentr[k] = ent 

               write( datout, string(conscore) ); write( datout, '\t' )
               write( datout, string(neutral) ); write( datout, '\t' )
               write( datout, string(slow) ); write( datout, '\t' )
               write( datout, string(deltag) ); write( datout, '\t' )
               write( datout, string(len) ); write( datout, '\t' )
               write( datout, string(ent) ); write( datout, '\t' )
               write( datout, string(meancov) ); write( datout, '\t' )
               write( datout, string(dist) ); write( datout, '\n' )
           
         end

         zcons = zscore(dcons)
         zslow = zscore(dslow)
         zdelg = zscore(ddelg)
         zdist = zscore(ddist)
         zentr = zscore(dentr)

         ztab  = DataFrame(Cons=zcons+zslow, DeltaG=zdelg, Dist=zdist)
         zmat  = transpose(Array(ztab))
         zdist = pairwise(Euclidean(), zmat)
         clust = dbscan(zdist, 1.5, 100)
         ztab[:,:cluster] = clust.assignments
         outliers = find(clust.assignments .== 0)

         for i in 1:length(clust.assignments)
            ind = indices[i]
            println(trav[ind])
            println("Z Conservation: $(ztab[i,:Cons])")
            println("Z DeltaG: $(ztab[i,:DeltaG])")
            println("Z Entropy: $(zentr[ind])")
            println("Delta G: $(ddelg[ind])")
            println("Conservation: $(dslow[ind])")
            println("Entropy: $(dentr[ind])")

            #if #ztab[i,:Cons] > 0.5 &&
               #ztab[i,:DeltaG] < 0 &&
            if dslow[ind] > 0.0 &&
               ddelg[ind] < -10 &&
               dentr[ind] > 0.12 

               println("OUTLIER!")
               writebed( bedout, trav[ind], "$(r.metadata):Outlier" )
            else
               ztab[i,:cluster] = 1
            end
         end

         # print plots
         consplot = plot(ztab, x=:DeltaG, y=:Cons, color=:cluster, Scale.color_discrete_manual("orange", "blue"), Theme(highlight_width = 0pt))
         distplot = plot(ztab, x=:DeltaG, y=:Dist, color=:cluster, Scale.color_discrete_manual("orange", "blue"), Theme(highlight_width = 0pt))
         draw(SVGJS("$(r.metadata)-Outlier.svg", 8inch, 3.5inch), hstack(consplot, distplot))
      end
   end
      
   close(datout)
   close(bedout)
end

@time main()


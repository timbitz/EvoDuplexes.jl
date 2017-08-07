#!/usr/bin/env julia

tic()
println( STDERR, "EvoDuplexes.jl $ver loading and compiling... " )

using ArgParse

unshift!( LOAD_PATH, dir * "/../src" )
using EvoDuplexes

function parse_cmd()
   s = ArgParseSettings()
   @add_arg_table s begin
      "tree"
       help     = "Phylogenetic tree with neutral branch lengths (in newick format)"
       arg_type = String
       required = true
      "regions"
       help     = "BED file containing regions for training"
       arg_type = String
       required = true
      "maf"
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

   regions = open(BEDReader, args["regions"])
   

end

main()


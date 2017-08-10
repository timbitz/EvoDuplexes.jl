#__precompile__()

module EvoDuplexes

using Bio
using Bio.Seq
using Bio.Intervals
using IntervalTrees
using SuffixArrays
using Gadfly
using Automa

using BufferedStreams
using Libz

import Automa
import Automa.RegExp: @re_str
import Compat: take!

importall Bio.Intervals

include("../src/pairs.jl")
include("../src/energy.jl")
include("../src/rnaduplex.jl")
include("../src/intervals.jl")
include("../src/traverse.jl")
include("../src/trie.jl")
include("../src/mafreader.jl")
include("../src/gtrmodel.jl")
include("../src/newick.jl")
include("../src/suffix.jl")
include("../src/evoduplex.jl")
include("../src/io.jl")
include("../src/regions.jl")

export RNADuplex,
       EvoDuplex,
       MAFSpecies,
       MAFRecord,
       MAFReader,
       readmaf!,
       read!,
       done,
       RNASuffixArray,
       RNADuplexArray,
       DuplexCollection,
       DuplexInterval,
       energy,
       npairs,
       traverse,
       collect,
       stitch,
       covariance,
       distance,
       score!,
       score,
       dinucleotide_entropy,
       writebed,
       extend_branches!,
       set_prob_mat!,
       parsenewick,
       GTR_SINGLE_Q,
       GTR_PAIRED_Q,
       loadbed


end

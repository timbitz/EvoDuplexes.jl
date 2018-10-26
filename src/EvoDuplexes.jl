#__precompile__()

module EvoDuplexes

using BioSymbols
using BioSequences
using GenomicFeatures
using IntervalTrees
using SuffixArrays
using Gadfly
using Automa

using BufferedStreams
using Libz

using PyCall
using PyPlot

using ScikitLearn
using ScikitLearn.Utils: meshgrid

@sk_import ensemble: IsolationForest

@pyimport matplotlib.font_manager as fm
@pyimport scipy.stats as stats

import Automa
import Automa.RegExp: @re_str
import Compat: take!

importall BioSymbols
importall GenomicFeatures

include("../src/pairs.jl")
include("../src/energy.jl")
include("../src/rnaduplex.jl")
include("../src/intervals.jl")
include("../src/traverse.jl")
include("../src/mafreader.jl")
include("../src/gtrmodel.jl")
include("../src/newick.jl")
include("../src/suffix.jl")
include("../src/evoduplex.jl")
include("../src/io.jl")
include("../src/regions.jl")
include("../src/train.jl")

export RNADuplex,
       EvoDuplex,
       MAFSpecies,
       MAFRecord,
       MAFReader,
       reverse_complement!,
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
       extend_branches!,
       set_prob_mat!,
       parsenewick,
       GTR_SINGLE_Q,
       GTR_PAIRED_Q,
       isgzipped,
       randcode,
       loadbed,
       writebed,
       loadbedgraph,
       DistanceForest,
       train!,
       predict
end

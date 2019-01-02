# EvoDuplexes.jl

[![Build Status](https://travis-ci.com/timbitz/EvoDuplexes.jl.svg?token=R7mZheNGhsReQ7hn2gdf&branch=master)](https://travis-ci.com/timbitz/EvoDuplexes.jl)
[![codecov](https://codecov.io/gh/timbitz/EvoDuplexes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/timbitz/EvoDuplexes.jl)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

RNA suffix array traversal through chromosomal multiple alignment files (MAF format) to fold and score all local and long-range RNA duplexes while allowing bulges and mismatches.

Full documentation is coming soon (**Requires Julia v0.6.4**), however comperehensive examples of API usage for EvoDuplexes package can be found in `test/runtests.jl`.

Executable can be run using `evo-duplex.jl` in /bin:
```bash
$ julia evo-duplex.jl -h
usage: evo-duplex.jl --tree TREE [--cons-regions CONS-REGIONS]
                     --gene-regions GENE-REGIONS --maf MAF
                     [--model-load MODEL-LOAD]
                     [--model-data MODEL-DATA] [--model-train]
                     [-o OUTPUT] [--cons-branch CONS-BRANCH]
                     [--sig-ratio SIG-RATIO] [--output-all]
                     [--max-distance MAX-DISTANCE]
                     [--max-duplex-length MAX-DUPLEX-LENGTH]
                     [--max-bulges MAX-BULGES]
                     [--max-mismatches MAX-MISMATCHES]
                     [--max-deltag MAX-DELTAG] [-h]

optional arguments:
  --tree TREE           Phylogenetic tree with neutral branch lengths
                        (in newick format)
  --cons-regions CONS-REGIONS
                        BED file containing conserved regions
  --gene-regions GENE-REGIONS
                        BED file containing gene units to allow
                        long-range folding within
  --maf MAF             Directory with MAF files named by chromosome
                        (chr1.maf.gz...) (default: "../maf")
  --model-load MODEL-LOAD
                        Load pre-trained IsolationForest models,
                        .evt.jls
  --model-data MODEL-DATA
                        Load training data from `.jlt` file, output
                        `.evt.jls` file
  --model-train         Train the IsolationForest models from scratch
                        using input regions, output `.evt.jls` file
  -o, --output OUTPUT   Output prefix (default: "evoduplex")
  --cons-branch CONS-BRANCH
                        When conserved regions are given, use this
                        branch length multiplier (type: Float64,
                        default: 0.33)
  --sig-ratio SIG-RATIO
                        Output the top fraction of outliers (type:
                        Float64, default: 0.05)
  --output-all          Output all folds regardless of prediction
                        status
  --max-distance MAX-DISTANCE
                        Set a limit on the maximum distance between a
                        left/right arm of a duplex (type: Int64,
                        default: 2000)
  --max-duplex-length MAX-DUPLEX-LENGTH
                        Set the maximum duplex length that will be
                        accessible from the suffix array (type: Int64,
                        default: 50)
  --max-bulges MAX-BULGES
                        Set the maximum number of bulges to allow for
                        a duplex (type: Int64, default: 3)
  --max-mismatches MAX-MISMATCHES
                        Set the maximum number of mismatches to allow
                        for a duplex (type: Int64, default: 3)
  --max-deltag MAX-DELTAG
                        Set the maximum allowable deltaG value (type:
                        Float64, default: -0.8)
  -h, --help            show this help message and exit
```

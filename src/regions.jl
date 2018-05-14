
function hasintersection( it::IntervalCollection{S}, i::Interval{T} ) where {S,T}
   !haskey( it.trees, i.seqname ) && return false
   try
      f = first(IntervalTrees.intersect(it.trees[i.seqname], i))
      return true
   catch
      return false
   end
end

function loadbed( filename::String; expandfirst::Int=0, expandlast::Int=0 )
   r  = BED.Reader(open( filename ))
   ic = IntervalCollection{String}()
   for b in r
      push!( ic, Interval{String}(seqname(b), 
                                  leftposition(b) - expandfirst, 
                                  rightposition(b) + expandlast,
                                  strand(b), String(b.data[b.name])) )
   end
   ic
end

# Eventually replace with Automa.jl parser...
function loadbedgraph( filename::String; regions::IntervalCollection{T}=IntervalCollection{Void}(), regionbool=false, expand=0 ) where T
   fh = open( filename , "r")
   if isgzipped( filename )
      fh = fh |> x->ZlibInflateInputStream(x, reset_on_end=true)
   end
 
   retcol = IntervalCollection{Float64}()

   for i in readlines(fh)
      ismatch( r"^#", i ) && continue
      spl = split(i, '\t')
      first = parse(Int, spl[2]) + 1
      last  = parse(Int, spl[3])
      if regionbool && !hasintersection( regions, Interval(spl[1], first-expand, last+expand) )
         continue
      else
         score = parse(Float64, spl[4])
         push!( retcol, Interval{Float64}(spl[1], first, last, '?', score) )
      end
   end

   retcol
end

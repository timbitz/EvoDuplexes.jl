
function hasintersection( it::IntervalCollection{S}, i::Interval{T} ) where {S,T}
   !haskey( it.trees, i.seqname ) && return false
   try
      f = first(IntervalTrees.intersect(it, i))
      return true
   catch
      return false
   end
end

function loadbed( filename::String; expandfirst::Int=0, expandlast::Int=0 )
   r  = BED.Reader(open(fixpath( filename )))
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
function loadbedgraph( filename::String; regions=IntervalCollection{T}, regionbool=false, expand=25 ) where T
   flat = fixpath( filename )
   fh = open( flat , "r")
   if isgzipped( flat )
      fh = fh |> x->ZlibInflateInputStream(x, reset_on_end=true)
   end
 
   retcol = IntervalCollection{Float64}()

   for i in readlines(hndl)
      ismatch( r"^#", i ) && continue
      spl = split(i, '\t')
      first = parse(Int, spl[2]) + 1
      last  = parse(Int, spl[3])
      if hasintersection( regions, Interval(spl[1], first-expand, last+expand) )
         score = parse(Float64, spl[4])
         push!( retcol, Interval{Float64}(spl[1], first, last, '?', score) )
      end
   end

   retcol
end

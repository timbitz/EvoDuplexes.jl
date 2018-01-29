
function loadbed( filename::String; expandfirst::Int=0, expandlast::Int=0 )
   r  = BED.Reader(open(filename))
   ic = IntervalCollection{String}()
   for b in r
      push!( ic, Interval{String}(seqname(b), 
                                  leftposition(b) - expandfirst, 
                                  rightposition(b) + expandlast,
                                  strand(b), String(b.data[b.name])) )
   end
   ic
end



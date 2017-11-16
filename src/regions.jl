
function loadbed( filename::String )
   r  = open(BEDReader, filename)
   ic = IntervalCollection{BEDMetadata}()
   for b in r
      push!( ic, b )
   end
   ic
end

function expandbed!( bed::IntervalCollection{BEDMetadata}, upstream::Int, downstream::Int )
   for r in bed
      r.first -= upstream
      r.last  += downstream
   end
   bed
end


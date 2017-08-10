
function loadbed( filename::String )
   r  = open(BEDReader, filename)
   ic = IntervalCollection{BEDMetadata}()
   for b in r
      push!( ic, b )
   end
   ic
end


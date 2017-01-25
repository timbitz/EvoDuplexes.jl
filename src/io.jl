
function tab_write( io, dat::Tuple )
   for i in 1:length(dat)
      write( io, string(dat[i]) )
      if i < length(dat)
         write( io, '\t' )
      end
   end
   write( io, '\n' )
end

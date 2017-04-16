
plug_write{S <: AbstractString}( io, str::S; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))
plug_write( io, str::Char; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))

tab_write{S <: AbstractString}( io, str::S ) = plug_write( io, str, plug='\t' )
tab_write( io, str::Char ) = plug_write( io, str, plug='\t' )

end_write{S <: AbstractString}( io, str::S ) = plug_write( io, str, plug='\n' )
end_write( io, str::Char ) = plug_write( io, str, plug='\n' )

function tab_write( io, dat::Tuple )
   for i in 1:length(dat)
      write( io, string(dat[i]) )
      if i < length(dat)
         write( io, '\t' )
      end
   end
   write( io, '\n' )
end

random_rgb() = (rand(0:255), rand(0:255), rand(0:255))
random_rgb_str() = join(random_rgb(), ",")

function writebed{T}( io, dup::DuplexInterval{T}, name::String; maxdistance::Int=100_000 )
   const dist = dup.last.first - dup.first.last
   if dup.first.seqname == dup.last.seqname &&
      dist <= maxdistance
      tab_write( io, string(dup.first.seqname) )
      tab_write( io, string(dup.first.first-1) )
      tab_write( io, string(dup.last.last) )
      tab_write( io, "DuplexTrie:" * name )
      tab_write( io, string(min(1000, energy(dup.duplex) * -250)) )
      tab_write( io, dup.first.strand == STRAND_POS ? '+' : '-' )
      tab_write( io, '.' )
      tab_write( io, '.' )
      tab_write( io, random_rgb_str() )
      tab_write( io, "2" )
      tab_write( io, string(length(dup.first.first:dup.first.last)) * "," * string(length(dup.last.first:dup.last.last)) )
      end_write( io, string(dup.first.first-1) * "," * string(dup.last.first-1) )
   end
end

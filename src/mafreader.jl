# Code based on Automa.jl multi-line 'FASTA' example:
# https://github.com/BioJulia/Automa.jl/blob/master/example/fasta.jl 

import Automa
import Automa.RegExp: @re_str
import Compat: take!
const re = Automa.RegExp

using Bio.Seq
using BufferedStreams


# Create a machine of MAF
const record_machine = (function ()

    lf          = re"\n"
    newline     = re"\r?" * lf
    lf.actions[:enter]          = [:countline]

    anything    = re.rep(re.diff( re".", newline ))
    header      = re.cat(re" ?#", anything, newline)

    float       = re.cat(re"[-+]?[0-9]*", re.opt("\."), re.opt(re"[0-9]*"))
    float.actions[:enter]       = [:mark]
    float.actions[:exit]        = [:score]

    sln         = re.rep1(re.space())
    emptyline   = re.cat(re" *", newline)

    scoreline   = re.cat(re" ?a score=", float, re.alt(newline, re.cat(re" ", anything, newline)))
    #scoreline.actions[:enter]   = [:setpos]

    organism    = re"[0-9A-Za-z\.\_\-]+"
    organism.actions[:enter]    = [:mark]
    organism.actions[:exit]     = [:organism]

    position    = re"[0-9]*"
    position.actions[:enter]    = [:mark]
    position.actions[:exit]     = [:position]

    size        = re"[0-9]*"
    size.actions[:enter]        = [:mark]
    size.actions[:exit]         = [:size]
 
    strand      = re"[+-]"
    strand.actions[:enter]      = [:mark]
    strand.actions[:exit]       = [:strand]

    sequence    = re"[ATGCNatgcn-]+"
    sequence.actions[:enter]    = [:mark]
    sequence.actions[:exit]     = [:sequence]

    entry       = re.cat(re" ?s", sln, organism, sln, position, sln, size, sln, strand, sln, re"[0-9]*", sln, sequence, newline)
    entry.actions[:enter]       = [:setpos]
    entry.actions[:exit]        = [:entry]

    fullheader  = re.cat(re.rep(header), re.rep(emptyline))
    fullheader.actions[:exit]   = [:anchor, :mark]

    record      = re.cat(re.opt(fullheader), re.opt(scoreline), entry, re.rep(entry), re.rep1(emptyline))
    record.actions[:enter]      = [:anchor]
    record.actions[:exit]       = [:escape]

    return Automa.compile(re.rep(record))
end)()

# It is useful to visualize the state machine for debugging.
# write("maf.dot", Automa.machine2dot(record_machine))
# run(`dot -Tsvg -o maf.svg maf.dot`)

# Bind Julia code to each action name (see the `parse_fasta` function defined below)

maf_actions_stream = Dict(
            :countline   => :(linenum += 1),
            :mark        => :(mark = p),
            :score => quote
                range = (mark:p-1) 
                score = length(range) <= 0 || mark == 0 ? "" : convert(String, data[range])
                mark = 0
            end,
            :organism => quote
                range = (mark:p-1)
                organism = mark == 0 ? "" : convert(String, data[range])
                mark = 0
            end,
            :position => quote
                range = (mark:p-1) 
                position = length(range) <= 0 || mark == 0 ? "" : convert(String, data[range])
                mark = 0
            end,
            :size => quote
                range = (mark:p-1) 
                size = length(range) <= 0 || mark == 0 ? "" : convert(String, data[range])
                mark = 0
            end,
            :strand => quote
                range = (mark:p-1) 
                strand = mark == 0 ? true : Char(data[last(range)]) == '+'
                mark = 0
            end,
            :sequence => quote
                range = (mark:p-1) 
                sequence = mark == 0 ? "" : chomp(String(data[range]))
                mark = 0
            end,
            :setpos => :(stream.position = p; mark_cs = cs), #println("Set stream position $p at cs $cs")),
            :entry  => quote
                try
                   #println(STDERR, "entry exit: $organism")
                   len    = parse(Int, size)
                   pscore = parse(Float64, score)
                   if len > reader.minlen && pscore > reader.minscore
                      push!(record, MAFSpecies(organism, parse(Int, position), strand, DNASequence(sequence)))
                   else
                      bad_record = true
                   end
                catch
                   println(STDERR, "size: $size and score: $score")
                   println(STDERR, "Bad Record at $p, setting p to $(stream.position)")
                   p = stream.position
                   @escape
                end
            end,
            :countline   => :(linenum += 1),
            :anchor      => :(stream.anchor = p), #; println("Anchoring at $p")),
            :escape      => :(found_record=true; @escape))

# Define a type to store a FASTA record.
type MAFSpecies
    name::String
    position::Int
    strand::Bool
    sequence::BioSequence{DNAAlphabet{4}}
end

Base.:(==)( a::MAFSpecies, b::MAFSpecies ) = a.name == b.name && a.position == b.position && 
                                             a.strand == b.strand && a.sequence == b.sequence ? true : false


type MAFReader{T<:BufferedInputStream}
    stream::T
    cs::Int
    p_eof::Int
    linenum::Int
    minlen::Int
    minscore::Float64
end

function MAFReader{BS <: BufferedInputStream}(input::BS)
    hits_eof = fillbuffer!(input) == 0
    p_eof = hits_eof ? input.available : -1
    return MAFReader{BS}(input, record_machine.start_state, p_eof, 1, 0, typemin(Float64))
end

Base.done(reader::MAFReader) = reader.cs <= 0

iscomment(line::String) = (length(line) >= 1 && line[1] == '#') || (length(line) >= 2 && line[2] == '#')
isheader(line::String)  = iscomment(line)

function discard_header!(reader::MAFReader)
    for l in eachline(reader.stream)
       reader.linenum += 1
       if isheader(l)
           continue
       end
       anchor!(reader.stream)
       hits_eof = fillbuffer!(reader.stream) == 0
       if hits_eof
          BufferedStreams.shiftdata!(reader.stream)
          reader.p_eof = reader.stream.available
       end
       reader.linenum -= 1
       return reader
    end
    reader.linenum -= 1
    reader
end

const MAFRecord = Vector{MAFSpecies}

@eval function Base.read!(reader::MAFReader, record::MAFRecord)
    cs      = reader.cs
    stream  = reader.stream
    data    = stream.buffer
    p       = stream.position
    p_end   = stream.available
    p_eof   = reader.p_eof
    linenum = reader.linenum
    mark = 0
    mark_cs = 0

    scoreline = organism = sequence = ""
    size = position = score = ""
    bad_record   = false
    found_record = false
    rebuffered   = false
    retry        = false

    if cs < 0
        error("the reader is in error state")
    elseif cs == 0
        error("the reader is finished")
    end

    while true

        p = p < 1 ? stream.anchor : p
        #error("p:$p p_end:$p_end p_eof:$p_eof cs:$cs anchor:$(stream.anchor) avail:$(stream.available) \"$(map(Char, data[1:10]))\"")

        $(Automa.generate_exec_code(record_machine, actions=maf_actions_stream, code=:goto, check=false))

        reader.cs              = cs
        reader.p_eof           = p_eof
        reader.linenum = linenum
      
        if found_record && !bad_record
           reader.stream.position = p
           return record
        elseif found_record && bad_record
           bad_record = false
           found_record = false
           continue
        elseif (p == p_eof+1 && data[p_eof] == 10)
           #println("p==p_eof+1 p:$p p_eof:$p_eof cs:$cs \"$(Char(data[p-1]))\" $record")
           spec = MAFSpecies(organism, parse(Int, position), strand, DNASequence(sequence))
           if !bad_record && !(record[end] == spec)
              push!(record, spec)
           end
           reader.cs = 0
           return record
        elseif p > p_eof ≥ 0
         #  println("p:$p p_eof:$p_eof cs:$cs \"$(Char(data[p-1]))\" $record")
           error("incomplete MAF input on line ", linenum)
        elseif p == p_end+1
            #println(STDERR, "Rebuffer p=p_end+1 p:$p p_eof:$p_eof p_end:$p_end cs:$cs")
            hits_eof = BufferedStreams.fillbuffer!(reader.stream) == 0
            if hits_eof
               reader.p_eof = p_eof = p_end
            else
               p     = stream.position
               p_end = stream.available
            end
            #println("Rebuffering to $(map(Char, data[p:p+10])) at cs $mark_cs")
            cs = record_machine.start_state
        elseif cs < 0
           if !retry && p <= 4
              retry = true
              p = 1
              cs = record_machine.start_state
           else
              println(String(data))
              println("cs:$cs p:$p p_end:$p_end $record anchor:$(stream.anchor) avail:$(stream.available) \"$(map(Char, data[1:p+10]))\"")
              error("MAF file format error on line ", linenum)
           end
        else
           error("Shouldn't happen cs:$cs p:$p p_eof:$p_eof $record")
        end

    end
end

function deletegaps!( mblock::MAFRecord )
   str = string(mblock[1].sequence)
   ind = search(str, r"-+")
   if length(ind) > 0
      deletegaps!( mblock, ind )
      deletegaps!( mblock )
      return
   else
      return
   end
end

function deletegaps!( mblock::MAFRecord, ind::UnitRange )
   for i in mblock
      deleteat!( i.sequence, ind )
   end
end


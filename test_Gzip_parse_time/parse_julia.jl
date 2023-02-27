using GZip
using TranscodingStreams
using CodecZlib

function parse_gzip_file(filename::String)
    out = open("test_gzipjl.txt", "w")
    zips = GZip.open(filename)
    while !eof(zips)
        line = readline(zips)
        println(out, line)
    end
    close(out)
end

function parse_gzi_trans(filename::String)
    out = open("test_trans.txt", "w")
    stream = GzipDecompressorStream(open(filename))
    for line in eachline(stream)
        println(out, line)
    end
    close(out)
end
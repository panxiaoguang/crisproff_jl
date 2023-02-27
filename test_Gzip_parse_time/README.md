### Compariation of time using for parsing same gzip file 

```
- python3                              0.399s
- julia -GZip.jl                       0.389s
- julia -TranscodingStreams.jl         0.280s
- shell                                0.297s
```
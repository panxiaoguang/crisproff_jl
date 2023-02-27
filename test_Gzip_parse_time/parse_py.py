import gzip
def parse_gzip_py(ris_file):
    inf = gzip.open(ris_file, 'rt')
    with open("test_py_gz.txt","w") as f:
        for line in inf:
            f.write(line+"\n")

    inf.close()
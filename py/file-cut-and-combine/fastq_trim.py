import os


def trim_umi(read, _max=11, _umi=6, _fix=3):
    """trim 6bp UMI and 3-5 Guanine"""
    start = _umi + _fix
    while start < _max:
        if read[start] != "G":
            break
        start += 1

    return start, read[:start], read[start:]


in_path = "/home/rcf-40/bos/3.scRNA-seq/download/fastq/"
# in_path = ""
with open(os.path.join(in_path, "filenames"), "r") as f:
    f.readline()
    fq_files = [line.rstrip() for line in f]

wrong_line_num = []
out_path = "/home/rcf-40/bos/3.scRNA-seq/download/trimmed_fastq"
# out_path = ""
# fq_files = ["test"]
for fq in fq_files:
    out = open(os.path.join(out_path, fq), "w")
    with open(os.path.join(in_path, fq), "r") as f:
        buff = ""
        for idx, line in enumerate(f):
            if idx % 4 == 1:
                trim_len, trimmed, reads = trim_umi(line)
                buff = (buff[:1] + trimmed + "_" + buff[1:].rstrip()
                        + "-" + str(trim_len) + "\n" + reads)
                continue
            if idx % 4 == 3:
                buff += line[trim_len:]
                out.write(buff)
                buff = ""
                continue
            buff += line
    out.close()

print(wrong_line_num)

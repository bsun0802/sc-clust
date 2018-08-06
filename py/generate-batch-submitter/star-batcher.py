from collections import defaultdict
from itertools import chain
import os

slurm = """#!/bin/sh
#SBATCH --partition=cmb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=51G
#SBATCH --time=200:00:00
#SBATCH --error=/home/rcf-40/bos/panfs/project/2.eQTL/log/Slurm/STAR_{}.e
#SBATCH --output=/home/rcf-40/bos/panfs/project/2.eQTL/log/Slurm/STAR_{}.o
#SBATCH --job-name=STAR_{}
"""
experiment = defaultdict(list)

with open("/home/rcf-40/bos/3.scRNA-seq/raw/run_info", "r") as f:
    f.readline()
    for line in f:
        rec = line.split()
        experiment[rec[8]].append(rec[6] + ".fastq")

fastq_path = "/home/rcf-40/bos/3.scRNA-seq/download/trimmed_fastq/"  # note the last /
c = 0
cmds = []
for exp in experiment:
    c += len(experiment[exp])
    readInFiles = ",".join([f"{fastq_path}{sra}" for sra in experiment[exp]])
    d = f"mkdir /home/rcf-40/bos/cmb00/3.scRNA-SEQ/SAM/{exp}"
    s = f"""/home/rcf-40/bos/panfs/software/STAR-2.6.0a/bin/Linux_x86_64/STAR \
--runThreadN 4 \
--outFilterScoreMinOverLread 0.55 \
--outFilterMatchNminOverLread 0.55 \
--outFileNamePrefix /home/rcf-40/bos/cmb00/3.scRNA-SEQ/SAM/{exp}/ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--genomeDir /home/rcf-40/bos/3.scRNA-seq/genome/M10_STAR_INDEX \
--readFilesIn {readInFiles}"""
    cmds.append((d, s))

target_sh_path = "/home/rcf-40/bos/panfs/project/2.eQTL/code/scRNA-seq/slurm/STAR"


def generate_sh(starter, cmds, i, path):
    out = open(os.path.join(path, f"{i // 450}_STAR_mapping_{i}"), "w")
    out.write(starter.format(i, i, i))
    out.write("\n".join(list(chain(*cmds))))
    out.close()


for i in range(0, 3000, 30):
    generate_sh(slurm, cmds[i:(i + 30)], i, target_sh_path)
    generate_sh(slurm, cmds[3000:3004], 3000, target_sh_path)

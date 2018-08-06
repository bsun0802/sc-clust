import os

gsm = []
with open(os.path.join("/Users/bos/lab/1.scRNA-seq/ncbi", "SraRunTable.txt"), "r") as f:
    f.readline()
    for line in f:
        rec = line.split()
        if rec[8] != "GSM1476865":
            gsm.append(rec[8])

slurm = """#!/bin/sh
#SBATCH --partition=cmb
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=200:00:00
#SBATCH --error=/home/rcf-40/bos/panfs/project/2.eQTL/log/Slurm/rsem_{}.e
#SBATCH --output=/home/rcf-40/bos/panfs/project/2.eQTL/log/Slurm/rsem_{}.o
#SBATCH --job-name=rsem_{}
"""


cmd = """rsem-calculate-expression --bam -p 8 \
/home/rcf-40/bos/staging/3.scRNA-seq/sam/{}/Aligned.toTranscriptome.out.bam \
/home/rcf-40/bos/cmb00/ref/M17PrimaryRSEM \
/home/rcf-40/bos/staging/3.scRNA-seq/sam/{}/rsem"""

cmds = list(set([cmd.format(id, id) for id in gsm]))

# print(cmds)
# print(len(cmds)) #3004


def generate_sh(starter, cmds, i, batch):
    out = open(os.path.join("/Users/bos/lab/1.scRNA-seq/code/batch-submitter-slurm",
                            f"{batch}_rsem_{i}"), "w")
    out.write(starter.format(i, i, i))
    out.write("\n".join(cmds))
    out.close()


batch = 1
for i in range(0, 2992, 34):
    generate_sh(slurm, cmds[i:(i + 34)], i, batch // 45)
    batch += 1

generate_sh(slurm, cmds[2992:3003], 2992, 1)

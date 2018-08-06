#!/usr/bin/env python3
# it's easy to go through the isoform expression file once to generate the gene and isoform id file
# in this case, os.path.join("/home/rcf-40/bos/3.scRNA-seq/output", "isoforms_id")
# contains the genes that have more than one isoform.
import os
import numpy as np

gsm = []  # 3004 cells
with open(os.path.join("/home/rcf-40/bos/3.scRNA-seq/raw", "run_info"), "r") as f:
    f.readline()
    for line in f:
        rec = line.split()
        if rec[8] != "GSM1476865":
            gsm.append(rec[8])

# all genes with >1 isoform will be kept
# gid = {}  # genes that have more than one isoform
# with open(os.path.join("/home/rcf-40/bos/3.scRNA-seq/output", "isoforms_id"), "r") as f:
#     f.readline()
#     for line in f:
#         gid[line.split()[0]] = 1


# cut out isoform TPM column and save them in _tmp
for acc_id in gsm:
    out = open(os.path.join(
        f"/home/rcf-40/bos/3.scRNA-seq/output/_tmp_all_tx_and_gene_tpm", f"{acc_id}_iso_TPM"), "w")
    out.write(acc_id + "\n")
    with open(os.path.join(f"/home/rcf-40/bos/staging/3.scRNA-seq/sam/{acc_id}",
                           "rsem.isoforms.results"), "r") as f:
        f.readline()
        for line in f:
            out.write(line.split()[5] + "\n")
    out.close()


# similarly, cut out all gene TPM and save them in _tmp
for acc_id in gsm:
    out = open(os.path.join(
        f"/home/rcf-40/bos/3.scRNA-seq/output/_tmp_all_tx_and_gene_tpm", f"{acc_id}_GENE_TPM"),
        "w")
    out.write(acc_id + "\n")
    with open(os.path.join(f"/home/rcf-40/bos/staging/3.scRNA-seq/sam/{acc_id}",
                           "rsem.genes.results"), "r") as f:
        f.readline()
        for line in f:
            out.write(line.split()[5] + "\n")
    out.close()


# Since isoformTPM.csv only contains the genes that have more than one transcript, So we will also need to cut the corresponding feature metadata.
# with open("/home/rcf-40/bos/staging/3.scRNA-seq/sam/GSM1476816/rsem.isoforms.results", "r") as f:
#     out = open("/home/rcf-40/bos/3.scRNA-seq/output/isoforms_meta.txt", "w")
#     out.write(f.readline())
#     for line in f:
#         rec = line.split()
#         if gid.get(rec[1]):
#             out.write(line)
#     out.close()


# with open("/home/rcf-40/bos/staging/3.scRNA-seq/sam/GSM1476816/rsem.genes.results", "r") as f:
#     out = open("/home/rcf-40/bos/3.scRNA-seq/output/genes_meta.txt", "w")
#     out.write(f.readline())
#     for line in f:
#         rec = line.split()
#         if gid.get(rec[0]):
#             out.write(line)
#     out.close()


# For helllinger distance, * deprecated *
# selecting completed cases only(i.e., isoform with 0 TPM across all cells will be ignored)
# cases = {}
# with open("/home/rcf-40/bos/3.scRNA-seq/output/tid.complete.cases", "r") as f:
#     f.readline()
#     while True:
#         line = f.readline()
#         if not line:
#             break
#         rec = line.split(",")
#         cases[(rec[0], rec[1])] = 1  # (tid, gid)


# def to_proportion(li):
#     a = np.array(li, dtype=float)
#     if np.sum(a):
#         a /= np.sum(a)
#     return [str(p) for p in a]


# for acc_id in gsm:
#     container = {}
#     with open(os.path.join(f"/home/rcf-40/bos/staging/3.scRNA-seq/sam/{acc_id}", "rsem.isoforms.results"), "r") as f:
#         f.readline()
#         for line in f:
#             rec = line.split()
#             if cases.get((rec[0], rec[1])):
#                 if container.get(rec[1]):
#                     container[rec[1]][0].append(rec[0])
#                     container[rec[1]][1].append(rec[5])
#                 else:
#                     container[rec[1]] = [[rec[0]], [rec[5]]]
#     out = open(f"/home/rcf-40/bos/3.scRNA-seq/output/_tmp_hellinger/{acc_id}_proportion", "w")
#     for k, v in container.items():
#         out.write(k + "\t" + ",".join(v[0]) + "\t" + ",".join(to_proportion(v[1])) + "\n")
#     out.close()

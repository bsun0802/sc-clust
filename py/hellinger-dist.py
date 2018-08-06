#!/usr/bin/env python3


# First attempt: reading in all isoform information and then filtering out genes that have more than one isoform as two dictionary, and then calculating the distance. With a one-time-writer idea.
# Deprecated reason: slow, each cell-file will be read in comparing every pair.
# def isoforms_expr(i, j) -> dict:
#     f1 = open("/Users/bos/lab/1.scRNA-seq/code/py/data/r2", "r")
#     f2 = open("/Users/bos/lab/1.scRNA-seq/code/py/data/rsem.isoforms.results", "r")
#     f1.readline()
#     f2.readline()
#     isoforms = {}
#     while True:
#         l1 = f1.readline()
#         l2 = f2.readline()
#         if not l1:
#             break
#         r1 = l1.split()
#         r2 = l2.split()
#         assert (r1[0] == r2[0] and r1[1] == r2[1])

#         if isoforms.get(r1[1]):
#             isoforms[r1[1]][0].append(r1[0])
#             isoforms[r1[1]][1].append(float(r1[5]))
#             isoforms[r1[1]][2].append(float(r2[5]))
#         else:
#             isoforms[r1[1]] = [[r1[0]], [float(r1[5])], [float(r2[5])]]

#     f1.close()
#     f2.close()
#     return isoforms


# def hellinger_distance(isoforms) -> dict:
#     global one_time_write
#     if one_time_write:
#         # for local test
#         out = open("/Users/bos/lab/1.scRNA-seq/code/py/output/isoforms_id", "w")
#         out.write("GENEID\tTRSPTID" + "\n")

#     isoforms_F = {}
#     for k, v in isoforms.items():
#         if len(v[0]) > 1:  # select out only genes that have more than 1 transcript
#             isoforms[k].append(to_proportion(v[1]))
#             isoforms[k].append(to_proportion(v[2]))
#             isoforms[k].append(h(isoforms[k][3], isoforms[k][4]))
#             isoforms_F[k] = isoforms[k]
#             if one_time_write:
#                 for tid in v[0]:
#                     out.write(f"{k}\t{tid}" + "\n")

#     if one_time_write:
#         out.close()
#         one_time_write = False
#     return isoforms_F
#

# ----


# Second Attempt: using a while loop to go through two files simultaneously.
# Deprecated reason: slower than attempt 1.
# NOTE: np.append is slower than list append, and then change to np.array
# def hellinger_distance(i, j) -> str:
#     H = 0
#     N = 0
#     # for local test
#     f1 = open("/Users/bos/lab/1.scRNA-seq/code/py/data/r2", "r")
#     f2 = open("/Users/bos/lab/1.scRNA-seq/code/py/data/rsem.isoforms.results", "r")
#     # f1 = open(os.path.join("/home/rcf-40/bos/staging/3.scRNA-seq/sam", i, "rsem.isoforms.results"), "r")
#     # f2 = open(os.path.join("/home/rcf-40/bos/staging/3.scRNA-seq/sam", j, "rsem.isoforms.results"), "r")
#     # Header = True
#     f1.readline()
#     f2.readline()
#     # init
#     rec1 = f1.readline().split()
#     rec2 = f2.readline().split()
#     curr_gene = rec1[1]
#     isoforms = [rec1[0]]
#     P = [float(rec1[5])]
#     Q = [float(rec2[5])]
#     while True:
#         l1 = f1.readline()
#         l2 = f2.readline()
#         if not l1:
#             break
#         r1 = l1.split()
#         r2 = l2.split()
#         if curr_gene == r1[1]:
#             isoforms.append(r1[0])
#             P.append(float(r1[5]))
#             Q.append(float(r2[5]))
#         else:
#             if len(isoforms) > 1:
#                 N += 1
#                 H += h(to_proportion(np.array(P)), to_proportion(np.array(Q)))
#             curr_gene = r1[1]
#             isoforms = [r1[0]]
#             P = [float(r1[5])]
#             Q = [float(r2[5])]

#     # Boundary case: the last gene need isn't taken care of by the while loop
#     if len(isoforms) > 1:
#         N += 1
#         H += h(to_proportion(np.array(P)), to_proportion(np.array(Q)))
#     f1.close()
#     f2.close()
#     print(N)
#     return str(H / N)

# ----


# Yet another deprecated approach:
# Read in all isoform proportion information into a large list of dictionaries, with one cell's isoforms proportion as a dict. Turned out it was too slow to use such a huge list.
# def isoform_tpm(acc_id) -> dict:
#     """{GENEID: [23.4, 0.0, 46.8...]}. TPM of each isoform"""
#     fn = os.path.join("/home/rcf-40/bos/staging/3.scRNA-seq/sam", acc_id, "rsem.isoforms.results")
#     ip = defaultdict(list)
#     with open(fn, "r") as f:
#         f.readline()
#         for line in f:
#             rec = line.split()
#             if gid.get(rec[1]):
#                 ip[rec[1]].append(float(rec[5]))
#     return ip


# gid = {}  # genes that have more than one isoform
# with open(os.path.join("/home/rcf-40/bos/3.scRNA-seq/output", "isoforms_id"), "r") as f:
#     f.readline()
#     for line in f:
#         gid[line.split()[0]] = 1

# ----


# Current approach: write isoform proportion of each cell to a file, and then build a huge list
import os
import numpy as np
# from collections import defaultdict
# from scipy.spatial.distance import euclidean


def catch_iso_proportion(acc_id) -> list:
    li_of_prop = []
    with open(os.path.join(_DATA, f"{acc_id}_proportion"), "r") as f:
        for line in f:
            li_of_prop.append(np.array(line.split()[-1].split(","), dtype=float))
    return li_of_prop


def _hellinger_distance(iso_i, iso_j) -> str:
    """Aggregated Hellinger's Distance between two cells, with each gene have the same weight"""
    H = 0
    for k in range(len(iso_i)):
        H += hellinger_distance(iso_i[k], iso_j[k])
    return str((H / _EQUAL_WEIGHT))


def hellinger_distance(P, Q) -> float:
    """Hellinger's Distance between two discrete probability distribution. Implemented as np.array of dtype=float.
    h(P, Q) = L2_norm( sqrt(P) - sqrt(Q) ) / sqrt(2).
    Divided by √2 is omitted for computational efficiency.

    NOTE: turned out that calculate the L-2 norm from scratch is faster than scipy.spatial.distance.euclidean, which is faster than np.linalg.norm
    """
    return np.sqrt(np.sum((np.sqrt(P) - np.sqrt(Q)) ** 2))


_EQUAL_WEIGHT = 19392  # 19392 genes, each with an equal weight
_LOCAL = False
_DATA = "/home/rcf-40/bos/3.scRNA-seq/output/_tmp_hellinger"
if _LOCAL:
    path = "/Users/bos/lab/1.scRNA-seq/code/py/data/"
    f1 = open(path + "r2_proportion", "r")
    f2 = open(path + "rsem.isoforms.results_proportion", "r")
    N = 0
    H = 0
    while True:
        l1 = f1.readline()
        l2 = f2.readline()
        if not l1:
            break
        N += 1
        H += hellinger_distance(np.array(l1.split()[-1].split(","), dtype=float),
                                np.array(l2.split()[-1].split(","), dtype=float))
    print(N, H / 20027 * (1 / np.sqrt(2)))
    f1.close()
    f2.close()
else:
    gsm = []  # accession ids
    with open(os.path.join("/home/rcf-40/bos/3.scRNA-seq/raw", "run_info"), "r") as f:
        f.readline()
        for line in f:
            rec = line.split()
            if rec[8] != "GSM1476865":
                gsm.append(rec[8])

    iso = [catch_iso_proportion(acc_id) for acc_id in gsm]

    with open(os.path.join("/home/rcf-40/bos/3.scRNA-seq/output/hellinger-distance", "dist-hellinger-new.txt"), "w") as distFile:
        prog = 0
        for i in range(len(gsm)):
            prog += 1
            for j in range(i + 1, len(gsm)):
                # write calculated Hellinger's distance of cell i and j to file
                distFile.write(
                    f"{gsm[i]}\t{gsm[j]}\t{_hellinger_distance(iso[i], iso[j])}" + "\n")
            print(f"The distance matrix has {prog} columns finished.")

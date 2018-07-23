# -*- coding: utf-8 -*-
"""
アミノ酸配列のFastaファイルを読み込んでペペペプチドでポジションを検索
アライメントの細かさを変えたい場合は--lengthや--mismatchの値を調整のこと
"""

from Bio import SeqIO
import requests
import argparse
import os.path

version = 3


def conPePePe(searchSeqs, mis):
    res = []

    pepepeUrl = "https://peptide.dbcls.jp/test/hg38_aa/"
    searchSeq = searchSeqs.pop(0)
    searchUrl = pepepeUrl + str(mis) + ":" + searchSeq
    for searchSeq in searchSeqs:
        searchUrl = searchUrl + " " + str(mis) + ":" + searchSeq
    searchUrl += ".txt"

    r = requests.get(searchUrl)
    rtxt = r.text
    rtxtarr = rtxt.split("\n")

    for line in rtxtarr:
        if line.startswith("#"):
            continue

        arr = line.split("\t")

        if len(arr) < 3:
            continue

        name = arr[0]
        length = int(arr[1])
        sta = int(arr[2])
        end = int(arr[3])
        mat = int(arr[7])
        mis = int(arr[8])

        chrid, frame = name.split("|")
        gsta = sta * 3
        gend = end * 3
        strand = "+"

        if frame == "ReadingFrame2":
            gsta += 1
        elif frame == "ReadingFrame3":
            gsta += 2
        elif frame == "ReadingFrame4":
            gsta = (length - end) * 3
            gend = (length - sta) * 3
            strand = "-"
        elif frame == "ReadingFrame5":
            gsta = (length - end) * 3 + 1
            gend = (length - sta) * 3 + 1
            strand = "-"
        elif frame == "ReadingFrame6":
            gsta = (length - end) * 3 + 2
            gend = (length - sta) * 3 + 2
            strand = "-"

        res.append([chrid, strand, gsta, gend, mat, mis])

    return res


def searchPePePe(infa, qlenDef, allowMis):
    outtxt, _ = os.path.splitext(infa)
    outtxt += "_pepepe.bed"
    outfh = open(outtxt, "w")

    for rec in SeqIO.parse(infa, "fasta"):
        aa = rec.seq
        aaid = rec.id
        aaseq = str(aa)
        aalen = len(aaseq)

        qind = 0
        chrstr2res = {}
        queries = []

        while qind < aalen:
            qlen = qlenDef

            if (aalen - qind - qlenDef * 2) < 0:
                qlen = aalen - qind
            elif (aalen - qind) < qlenDef:
                qlen = aalen - qind

            qseq = aaseq[qind:(qind+qlen)]
            queries.append(qseq)
            qind += qlen

        res = conPePePe(queries, allowMis)

        for r in res:
            chrid, strand, gsta, gend, mat, _ = r
            chrstr = chrid + strand
            if chrstr not in chrstr2res:
                chrstr2res[chrstr] = []

            name = aaid
            r.append(name)
            chrstr2res[chrstr].append(r)

        for chrstr, res in chrstr2res.items():
            for r in res:
                chrid, strand, gsta, gend, mat, _, name = r

                # bed of each alignment
                bed = [chrid, str(gsta - 1), str(gend), name, str(mat), strand]
                outfh.write("\t".join(bed) + "\n")

    outfh.close()


def main():
    par = argparse.ArgumentParser(description="align AA fasta on pepepeptide")

    par.add_argument("-i", "--input", required=True,
                     help="a fasta file of amino acids")

    par.add_argument("-l", "--length", type=int, default=10,
                     help="default length of query")

    par.add_argument("-m", "--mismatch", type=int, default=2,
                     help="number of mismatches")

    args = par.parse_args()
    infa = args.input
    qlenDef = args.length
    allowMis = args.mismatch

    searchPePePe(infa, qlenDef, allowMis)


if __name__ == '__main__':
    main()

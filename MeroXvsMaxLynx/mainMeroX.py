import pandas as pd
import MeroxID as merox
from Bio import SeqIO

import sys
sys.path.append('/Figures/MeroXvsMaxLynx')


### Collect protein sequences and their accession numbers
## First get accession and sequence information from ID.fasta
fasta = "C:/yilmaz/Revision/param_scanning/large_scale/MeroX_paper/ID.fasta"
idAseq = dict()
for seq_record in SeqIO.parse(fasta, "fasta"):
    identifier = seq_record.id
    acc = identifier.split("|")[1]
    idAseq[acc] = seq_record.seq

### Now collect PSMs on MeroX list and convert these to Link objects
meroxFile = "C:/yilmaz/Revision/param_scanning/large_scale/MeroX_paper/MeroX_toreport/MeroX_crosslinktable.csv"
df = pd.read_csv(meroxFile, sep=';')
allMeroxPsm = list()
meroX_Links = list()
meroX_IntraLinks = list()

print("Started.")
for index, row in df.iterrows():
    meroxPsm = merox.MeroxID(row.Score, row[1], row.Charge, row[3], row[4], row[5], row[6], row[7], row[8], row[9],
                             row[10], row[11], row[12], row[13], row[14], row[18], row[20], row[21],
                             row.pepScore1, row.pepScore2, row.xLinkScore, idAseq)
    allMeroxPsm.append(meroxPsm)
    if not meroxPsm.link in meroX_Links:
        meroX_Links.append(meroxPsm.link)

    if meroxPsm.link.xlsites1.keys() == meroxPsm.link.xlsites2.keys() and not meroxPsm.link in meroX_IntraLinks:
        meroX_IntraLinks.append(meroxPsm.link)
    # keyCount = 0
    #
    # for key1 in meroxPsm.link.xlsites1.keys():
    #     for key2 in meroxPsm.link.xlsites2.keys():
    #         if key1 == key2:
    #             keyCount = keyCount + 1
    # if keyCount == meroxPsm.link.xlsites1.keys():


print("Total MeroX_Links= ", str(len(meroX_Links)))
print("Total MeroX_IntraLinks= ", str(len(meroX_IntraLinks)))

with open('C:\yilmaz\Revision\param_scanning\large_scale\pcharm_output\meroX_SYR_Links_PyCharm_2012.txt', 'w') as f:
    f.write("Pro1;Pro2;Link1;Link2 \n")
    for t in meroX_Links:
        f.write(str(list(t.xlsites1.keys())))
        f.write(";")
        f.write(str(list(t.xlsites1.values())))
        f.write(";")
        f.write(str(list(t.xlsites2.keys())))
        f.write(";")
        f.write(str(list(t.xlsites2.values())))
        f.write('\n')
    f.close()

with open('C:\yilmaz\Revision\param_scanning\large_scale\pcharm_output\meroX_SYR_MeroXPSMs_PyCharm_2012.txt', 'w') as f:
    f.write(
        "Score;MZ;Peptide1;Proteins1;Peptide2;Proteins2;From1;From2;Site1;Site2;MeroXScan;LinkedPro1;LinkedProSites1;LinkedPro2;LinkedProSites2;Replicate \n")
    for t in allMeroxPsm:
        f.write(str(t.score))
        f.write(";")
        f.write(str(t.mz))
        f.write(";")
        f.write(t.peptide1)
        f.write(";")
        f.write(str(t.proteins1))
        f.write(";")
        f.write(t.peptide2)
        f.write(";")
        f.write(str(t.proteins2))
        f.write(";")
        f.write(str(t.from1))
        f.write(";")
        f.write(str(t.from2))
        f.write(";")
        f.write(str(t.bestLinkagePeptide1))
        f.write(";")
        f.write(str(t.bestLinkagePeptide2))
        f.write(";")
        f.write(t.meroxScan)
        f.write(";")
        f.write(str(list(t.link.xlsites1.keys())))
        f.write(";")
        f.write(str(list(t.link.xlsites1.values())))
        f.write(";")
        f.write(str(list(t.link.xlsites2.keys())))
        f.write(";")
        f.write(str(list(t.link.xlsites2.values())))
        f.write(";")
        f.write(t.replicate)
        f.write('\n')
    f.close()

import pandas as pd
import MaxLynxID as maxlynx
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

### Now collect PSMs on MaxLynx list and convert these to Link objects
settingsName = "0_10_SepFDR_FDR100"
maxlynxFile = "C:/yilmaz/Revision/param_scanning/large_scale/MaxLynx_toreport/crosslinkMsms_"+ settingsName + ".txt"
maxlynxFile = "C:/yilmaz/Revision/param_scanning/large_scale/fdr100/crosslinkMsms_0_10_SepFDR_FDR100.txt"
pycharmOutputEnd = settingsName + "_pycharm.txt"
print(maxlynxFile)

df = pd.read_csv(maxlynxFile, sep='\t')
allMaxLynxPsms = list()
maxLynxLinks = list()
intraLinks = list ()

print("Started.")
for index, row in df.iterrows():
    decoyType = row.Decoy
    if decoyType == "forward" and (row[13] == "IntraProXL" or row[13] == "InterProXL"):
        maxlynxPsm = maxlynx.MaxLynxID(row[0], row[1], row.Sequence1, row.Sequence2, row[5], row[6],
                 decoyType, row.Score, row.PEP, row[13], row.Proteins1, row.Proteins2, row[20], row[23], row[21], row[24],
                 row[43], row[44], row[45], row.Charge, row[47], row.Mass, row[49], row[50], idAseq)
        allMaxLynxPsms.append(maxlynxPsm)
        if not maxlynxPsm.link in maxLynxLinks:
            maxLynxLinks.append(maxlynxPsm.link)
        if row[13] == "IntraProXL" and not maxlynxPsm.link in intraLinks:
            intraLinks.append(maxlynxPsm.link)

print("#CSMS=", str(len(allMaxLynxPsms)))
print("#MaxLynxLinks= ", str(len(maxLynxLinks)))
print("#MaxLynx intraLinks= ", str(len(intraLinks)))
print("Ended.")


with open('C:\yilmaz\Revision\param_scanning\large_scale\pcharm_output\maxlynx_Links_' + pycharmOutputEnd, 'w') as f:
    f.write("Pro1;Pro2;Link1;Link2 \n")
    for t in maxLynxLinks:
        f.write(str(list(t.xlsites1.keys())))
        f.write(";")
        f.write(str(list(t.xlsites1.values())))
        f.write(";")
        f.write(str(list(t.xlsites2.keys())))
        f.write(";")
        f.write(str(list(t.xlsites2.values())))
        f.write('\n')
    f.close()


with open('C:\yilmaz\Revision\param_scanning\large_scale\pcharm_output\maxlynx_PSMs_' + pycharmOutputEnd, 'w') as f:
    f.write(
        "Score;MZ;Peptide1;Proteins1;Peptide2;Proteins2;LinkPosition1;LinkedPosition2;Site1;Site2;Scan;RawFile;LinkedPro1;LinkedProSites1;LinkedPro2;LinkedProSites2;Replicate\n")
    for t in allMaxLynxPsms:
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
        f.write(str(t.pepInterLink1))
        f.write(";")
        f.write(str(t.pepInterLink2))
        f.write(";")
        f.write(str(t.aa1))
        f.write(";")
        f.write(str(t.aa2))
        f.write(";")
        f.write(str(t.scanNum))
        f.write(";")
        f.write(t.rawFile)
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
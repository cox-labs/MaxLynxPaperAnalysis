## This class collects all the attributes of MeroX to compare against MQ output later
## Cross-linking on protein N-term is 1!! Now, also with MeroX result
## Be careful when you match proteins between MeroX and MaxLynx: Both protein-list are sorted within their lists, with key=str.lower
## Here all MaxLynx di-peptide identifications have only one linked residue per peptide.
## A linked residue with bigger value is written first.

import sys
sys.path.append('/Figures/MeroXvsMaxLynx')

import Link as ld


class MaxLynxID:


    ## This version of MaxLynx (MaxLynx2.0.2-2.0.3), all proteins were reported within paranthesis.
    ## Here any decoys were excluded (similar to MeroX output)
    def GetProteins(self, proteins):
        maxLynxProteins = list()
        if "(" in proteins:
            prosArr = proteins.split("(")[1].split(")")[0].split(";")
            for pros in prosArr:
                if not pros.startswith("REV__") and not "|" in pros:
                    maxLynxProteins.append(pros)
                if not pros.startswith("REV__") and  "|" in pros:
                    maxLynxProteins.append(pros.split("|")[1])
        else:
            if not "|" in proteins:
                maxLynxProteins.append(proteins)
            else:
                maxLynxProteins.append(proteins.split("|")[1])
        maxLynxProteins = sorted(maxLynxProteins, key=str.lower)
        return maxLynxProteins

    ## From given proteins, peptide, idAseq and site information,
    ## this method finds which proteins contains a given peptide and based on siteInfo it gives a linking position on this protein
    ## This returns a dictionary: protein-names as a key and linking-position on this protein as a value
    ## Protein n-term cross-linking (meaning any crosslinked residues except K) has the linking site of 1 (similar to MeroX)
    def GetXLSites(self, proteins, peptide, idAseq, linkPositionPeptide, linkResidue):
        proAind = {}
        for tmpProtein in proteins:
            index = idAseq[tmpProtein].find(peptide)
            if "K" in linkResidue:
                siteInt = int(linkPositionPeptide)
                tmpIndex = index + siteInt
            else:
                tmpIndex = 1
            proAind[tmpProtein] = tmpIndex
        return proAind

    def GetLinkedIndexForFirstProtein(self, proteins, peptide, idAseq, linkPosition, linkResidue):
        startIndex = idAseq[proteins[0]].find(peptide)
        if "K" in linkResidue:
            siteInt = int(linkPosition)
            linkedIndex = startIndex + siteInt
        else:
            linkedIndex = 1
        return linkedIndex



    def __init__(self, rawFile, scanNum, peptide1, peptide2, modPeptide1, modPeptide2,
                 decoy, score, pep, product, proteins1, proteins2, pepInterLink1, pepInterLink2, aa1, aa2,
                 deltaScore, partialScore1, partialScore2,
                 charge, mz, mass, ppmerr, annotations, idAseq):

        self.rawFile, self.scanNum = rawFile, scanNum
        self.replicate = rawFile.split("_")[0]
        self.peptide1, self.peptide2, self.modPeptide1, self.modPeptide2 = peptide1, peptide2, modPeptide1, modPeptide2
        self.aa1, self.aa2, self.pepInterLink1, self.pepInterLink2 = aa1, aa2, pepInterLink1, pepInterLink2
        if ";" in aa1:
            self.aa1 = aa1[:-1]
            self.pepInterLink1 = pepInterLink1[:-1]
        if ";" in aa2:
            self.aa2 = aa2.split(";")[0]
            self.pepInterLink2 = pepInterLink2[:-1]
        self.proteins1, self.proteins2 = self.GetProteins(proteins1), self.GetProteins(proteins2)

        self.decoy = decoy
        self.score, self.pep = score, pep
        self.product = product

        self.deltaScore = deltaScore
        self.partialScore1, self.partialScore2 = partialScore1, partialScore2
        self.charge, self.mz, self.mass, self.ppmerr = charge, mz, mass, ppmerr
        self.annotations = annotations

        self.xlSites1 = self.GetXLSites(self.proteins1, self.peptide1, idAseq, self.pepInterLink1, self.aa1)
        self.xlSites2 = self.GetXLSites(self.proteins2, self.peptide2, idAseq, self.pepInterLink2, self.aa2)

        self.link = ld.Link(self.xlSites1, self.xlSites2)

        ## To keep the same order to reporting crosslinks as well as avoiding redundant cross-links, check the link-position of the first element on each link-list and report a crosslink-pair with smaller index-value as the first-link-site (as what MeroX does as well)
        linkSites1FirstProtein = self.GetLinkedIndexForFirstProtein(self.proteins1, self.peptide1, idAseq, self.pepInterLink1, self.aa1)
        linkSites2FirstProtein = self.GetLinkedIndexForFirstProtein(self.proteins2, self.peptide2, idAseq, self.pepInterLink2, self.aa2)

        if int(linkSites1FirstProtein) > int(linkSites2FirstProtein):
            self.link = ld.Link(self.xlSites2, self.xlSites1)

    def __str__(self):
        return ("rawFile: ", self.rawFile, " scanNum: ", str(self.scanNum), " replicate: ", str(self.replicate),
                " Peptide1: ", self.peptide1, " Peptide2: ", self.peptide2, " ModPeptide1: ", self.modPeptide1,
                " ModPeptide2: ", self.modPeptide2,
                " Decoy: ", self.decoy, " Score: ", str(self.score), " PEP: ", str(self.pep),
                " ProductType: ", self.product, " Proteins1: ", self.proteins1, " Proteins2: ", self.proteins2,
                " Link1: ", str(self.pepInterLink1), " Link2: ", str(self.pepInterLink2),
                " aa1: ", self.aa1, " aa2: ", self.aa2, " DeltaScore: ", str(self.deltaScore), " PartialScore1: ",
                str(self.partialScore1), " PartialScore2: ", str(self.partialScore2),
                " Charge: ", str(self.charge), " MZ: ", str(self.mz), "Mass: ", str(self.mass), " ppm-error: ",
                str(self.ppmerr), " Annotations: ", self.annotation,
                " Link: ", self.link.toString())

## This class collects all the attributes of MeroX (which was exported to csv file from the MeroX session)
## Cross-linking event on protein N-term is 0 for the original MeroX but now it is 1, similar to MaxLynx!!
## Be careful when you match proteins between MeroX and MaxLynx: Proteins are sorted within their lists, with key=str.lower

import sys
sys.path.append('/Figures/MeroXvsMaxLynx')

import Link as ld


class MeroxID:

    ## MeroX pre-processed already protein description in case that there is "(".
    ## For example, ID.fasta file > ">sp|Q9V3N8|THG1_DROME Probable tRNA(His) guanylyltransferase OS=Drosophila melanogaster OX=7227 GN=THG PE=1 SV=1" became as follows:
    ## ">sp|Q9V3N8|THG1_DROME Probable tRNA His  guanylyltransferase OS=Drosophila melanogaster OX=7227 GN=THG PE=1 SV=1"
    def GetFirstMeroXProtein(self, proteins):
        if "(" in proteins:
            prosArr = proteins.split("(")
            firstProtein = prosArr[0].split("|")[1]
            return firstProtein
        firstProtein = proteins.split("|")[1]
        return firstProtein

    ## A MeroX protein-list example: >sp|Q99323|MYSN_DROME Myosin heavy chain, non-muscle OS=Drosophila melanogaster OX=7227 GN=zip PE=1 SV=2(>tr|A0A0B4JD57|A0A0B4JD57_DROME Zipper, isoform F OS=Drosophila melanogaster OX=7227 GN=zip PE=1 SV=1/>tr|A0A0B4JD95|A0A0B4JD95
    ## Only protein accession-numbers (after being sorted) are returned.
    def GetProteins(self, proteins):
        proteinList = list()
        if "(" in proteins:
            prosArr = proteins.split("(")
            restPros = prosArr[1][:-1].split("/")
            for tmp in restPros:
                proteinList.append(tmp.split("|")[1])
            proteinList.append(prosArr[0].split("|")[1])
        else:
            proteinList.append(proteins.split("|")[1])
        proteinList = sorted(proteinList, key=str.lower)
        return proteinList

    ## From given proteins, peptide, idAseq and site information,
    ## this method finds which proteins contains a given peptide and based on siteInfo it gives a linking position on this protein
    ## This returns a dictionary: protein-names as a key and linking-position on this protein as a value
    ## Protein n-term cross-linking has the linking site of 1
    def GetXLSites(self, proteins, peptide, idAseq, site):
        proAind = {}
        for tmpProtein in proteins:
            index = idAseq[tmpProtein].find(peptide)
            if "K" in site:
                siteInt = int(site.split("K")[1])
                tmpIndex = index + siteInt
            elif "{0":
                tmpIndex = 1
            proAind[tmpProtein] = tmpIndex
        return proAind

    def GetLinkedIndexForFirstProtein(self, proteins, peptide, idAseq, site):
        startIndex = idAseq[proteins[0]].find(peptide)
        if "K" in site:
            siteInt = int(site.split("K")[1])
            linkedIndex = startIndex + siteInt
        elif "{0":
            linkedIndex = 1
        return linkedIndex

    ## Depending on the MeroX scan information, this method returns the isotope shift-value done by MeroX as well as spectrumName (i.e. raw file name)
    def GetIsotopeShiftAspecName(self, scan):
        specInfo = scan.split(" ")
        if len(specInfo) == 2:
            isotopeShift = int(specInfo[0][1:-1])
            specName = specInfo[1].split("~")
        else:
            isotopeShift = 0
            specName = specInfo[0].split("~")
        return isotopeShift, specName

    def __init__(self, score, mz, charge, precMass, calcMass, massErr,
                 peptide1, proteins1, from1, to1,
                 peptide2, proteins2, from2, to2,
                 scan, rt,
                 bestLinkagePeptide1, bestLinkagePeptide2,
                 pepScore1, pepScore2, xLinkScore,
                 idAseq):

        self.score, self.mz, self.charge, self.precMass, self.calcMass, self.massErr = score, mz, charge, precMass, calcMass, massErr

        self.meroxPeptide1, self.meroxProteins1, self.from1, self.bestLinkagePeptide1 = peptide1, proteins1, from1, bestLinkagePeptide1
        # Assuring to find the same peptide on the given fasta file: MeroX peptides are written within [], oxidation is shown as "m" (GKmQPTHPIR ==> m is oxidized) and C>B (as fixed/static modification)
        self.peptide1 = peptide1[1:-1].replace('B', 'C').upper()
        self.firstProtein1, self.proteins1 = self.GetFirstMeroXProtein(proteins1), self.GetProteins(proteins1)
        self.xlSites1 = self.GetXLSites(self.proteins1, self.peptide1, idAseq, bestLinkagePeptide1)

        self.meroxPeptide2, self.meroxProteins2, self.from2, self.bestLinkagePeptide2 = peptide2, proteins2, from2, bestLinkagePeptide2
        self.peptide2 = peptide2[1:-1].replace('B', 'C').upper()
        self.firstProtein2, self.proteins2 = self.GetFirstMeroXProtein(proteins2), self.GetProteins(proteins2)
        self.xlSites2 = self.GetXLSites(self.proteins2, self.peptide2, idAseq, bestLinkagePeptide2)

        self.link = ld.Link(self.xlSites1, self.xlSites2)

        ## To keep the same order to reporting crosslinks as well as avoiding redundant cross-links, check the link-position of the first element on each link-list and report a crosslink-pair with smaller index-value as the first-link-site (as how MeroX-original reports as well)
        linkSites1FirstProtein = self.GetLinkedIndexForFirstProtein(self.proteins1, self.peptide1, idAseq, bestLinkagePeptide1)
        linkSites2FirstProtein = self.GetLinkedIndexForFirstProtein(self.proteins2, self.peptide2, idAseq, bestLinkagePeptide2)
        if int(linkSites1FirstProtein) > int(linkSites2FirstProtein):
            self.link = ld.Link(self.xlSites2, self.xlSites1)

        self.meroxScan, self.rt = scan, rt
        self.isotopeShift, specName = self.GetIsotopeShiftAspecName(scan)
        self.scanNum, self.rawFile = specName[0], specName[1]
        self.replicate = specName[1].split("_")[0]

        self.pepScore1, self.pepScore2, self.xLinkScore = pepScore1, pepScore2, xLinkScore

    def __str__(self):
        return ("score: ", str(self.score), " pepscore1: ", str(self.pepScore1), " pepscore2: ", str(self.pepScore2),
                " mz: ", str(self.mz), " Charge: ", str(self.charge), " calc: ", str(self.calcMass), " mass-error:",
                str(self.massErr),  " MeroXPeptide1: ", self.meroxPeptide1, " peptide1: ", self.peptide1, " meroxproteins1: ",
                self.meroxProteins1, " proteins1: ", self.proteins1,  " firstProtein1: ", self.firstProtein1, " xlsites1: ", str(self.xlSites1),
                " MeroXPeptide2: ", self.meroxPeptide2, " peptide2: ", self.peptide2, " meroxproteins2: ",
                self.meroxProteins2, " proteins2: ", self.proteins2, " firstProtein2: ", self.firstProtein2, " xlsites2: ", str(self.xlSites2),
                " Scan: ", str(self.meroxScan),  " RawFile: ", self.rawFile, " Scannum: ", str(self.scanNum), " Replicate: ", self.replicate,
                " IsotopeShift: ", str(self.isotopeShift), " RT: ", str(self.rt), " Link: ", str(self.link))

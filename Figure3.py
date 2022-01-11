## PLOT FDR=1% CSMs level - Note that original data did not provide Kojak, so it only available for FDR=5% data
import pandas as pandas
import matplotlib as matplotlib
import matplotlib.pyplot as plt

df = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'],  correct=[737, 940, 880], incorrect=[6, 18, 11]))
dfOpen = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[368, 506, 365], incorrect=[4, 5, 4]))
dfpLink = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[594, 644, 585], incorrect=[10, 13, 25]))
dfStavroX = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[265, 157, 160], incorrect=[4, 0, 1]))
dfXi = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[312, 352, 438], incorrect=[2, 4, 5]))

matplotlib.rc('ytick', labelsize=10)
matplotlib.rc('xtick', labelsize=10)

def plotForManuscript(pdDfTotal, pdfCorrect, increase, programName, programNameLine, writeLegend):
    ind = [0.5 + increase, 1 + increase, 1.5 + increase]
    width = 0.5
    ax.barh(ind, pdDfTotal, width, color='#e74c3c', label='incorrect', alpha=1, edgecolor='white', linewidth=0.4)
    ax.barh(ind, pdfCorrect, width, color='#86c5da', label='correct', alpha=1, edgecolor='white', linewidth=0.4)
    ax.text(programNameLine, 0.75 + increase, programName, verticalalignment='bottom', rotation=0, multialignment='right')
    if (writeLegend):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right')

fig, ax = plt.subplots()
ax.set_title("FDR=1%, #CSMs")
plotForManuscript(df.incorrect + df.correct, df.correct, 0, "MaxLynx", -300, True)
plotForManuscript(dfOpen.incorrect + dfOpen.correct, dfOpen.correct, 2, "OpenPepXL", -350, False)
plotForManuscript(dfpLink.incorrect + dfpLink.correct, dfpLink.correct, 4, "pLink2", -250, False)
plotForManuscript(dfStavroX.incorrect + dfStavroX.correct, dfStavroX.correct, 6, "StavroX", -270, False)
plotForManuscript(dfXi.incorrect + dfXi.correct, dfXi.correct, 8, "Xi", -150, False)

ylabels = df.graph.append(df.graph).append(df.graph).append(df.graph).append(df.graph)
inds = [0.5, 1, 1.5, 2.5, 3, 3.5, 4.5, 5, 5.5, 6.5, 7, 7.5, 8.5, 9, 9.5]
ax.set(yticks=inds, yticklabels=ylabels, xlim=[0, 1500], ylim=[10, 0])
csmfdr001 = ax

plt.savefig('C:/yilmaz/Revision/Figures/Figure3a_dss_csms_fdr001..png', dpi=300, bbox_inches='tight')

## PLOT FDR=1% crosslink level - Note that original data did not provide Kojak, so it only available for FDR=5% data

df = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[227, 240, 223], incorrect=[6, 14, 9]))
dfOpen = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[161, 196, 148], incorrect=[2, 4, 3]))
dfpLink = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[215, 218, 189], incorrect=[9, 12, 22]))
dfStavroX = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[124, 91, 90], incorrect=[4, 0, 1]))
dfXi = pandas.DataFrame(dict(graph=['R1', 'R2', 'R3'], correct=[141, 152, 163], incorrect=[2, 3, 5]))
matplotlib.rc('ytick', labelsize=10)
matplotlib.rc('xtick', labelsize=10)

def plotForManuscript(pdDfTotal, pdfCorrect, increase, programName, programNameLine, writeLegend):
    ind = [0.5 + increase, 1 + increase, 1.5 + increase]
    width = 0.5
    ax.barh(ind, pdDfTotal, width, color='#e74c3c', label='incorrect', linewidth=0.4, alpha=1, edgecolor='white')
    ax.barh(ind, pdfCorrect, width, color='#2980b9', label='correct', linewidth=0.4, alpha=1, edgecolor='white')
    ax.text(programNameLine, 0.75 + increase, programName, verticalalignment='bottom', rotation=0, multialignment='right')
    if (writeLegend):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right')

fig, ax = plt.subplots()
ax.set_title("FDR=1%, #crosslinks")
plotForManuscript(df.incorrect + df.correct, df.correct, 0, "MaxLynx", -80, True)
plotForManuscript(dfOpen.incorrect + dfOpen.correct, dfOpen.correct, 2, "OpenPepXL", -95, False)
plotForManuscript(dfpLink.incorrect + dfpLink.correct, dfpLink.correct, 4, "pLink2", -66, False)
plotForManuscript(dfStavroX.incorrect + dfStavroX.correct, dfStavroX.correct, 6, "StavroX", -75, False)
plotForManuscript(dfXi.incorrect + dfXi.correct, dfXi.correct, 8, "Xi", -40, False)

ylabels = df.graph.append(df.graph).append(df.graph).append(df.graph).append(df.graph)
inds = [0.5, 1, 1.5, 2.5, 3, 3.5, 4.5, 5, 5.5, 6.5, 7, 7.5, 8.5, 9, 9.5]
ax.set(yticks=inds, yticklabels=ylabels, xlim=[0, 400], ylim=[10, 0])
crosslinkfdr001 = ax

plt.savefig('C:/yilmaz/Revision/Figures/Figure3b_dss_crosslink_fdr001.png', dpi=300, bbox_inches='tight')
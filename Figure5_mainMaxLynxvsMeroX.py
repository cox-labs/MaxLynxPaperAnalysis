## This class reads the output from mainmerox.py and mainmaxlynx.py classes
import sys
sys.path.append('/MeroXvsMaxLynx')

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2
import seaborn as sns

### Put these three venn diagrams next to each other

### Now collect PSMs on MeroX list and convert these to Link objects
############ Reading all MeroX data ###########
print("Started:  reading MeroX and MaxLynx outputs")
meroxFile = "C:/yilmaz/Revision/param_scanning/large_scale/pcharm_output/meroX_SYR_MeroXPSMs_PyCharm_2012.txt"
df = pd.read_csv(meroxFile, sep=';')
meroX_Links, meroXLinksReplicate1, meroXLinksReplicate2, meroXLinksReplicate3 =[],[],[],[]

for index, row in df.iterrows():
    tmpReplicate = row.Replicate
    tmpLink = row.LinkedPro1, row.LinkedProSites1, row.LinkedPro2, row.LinkedProSites2
    if not tmpLink in meroX_Links:
        meroX_Links.append(tmpLink)
    if row.Replicate == "R1" and not tmpLink in meroXLinksReplicate1:
        meroXLinksReplicate1.append(tmpLink)
    elif row.Replicate == "R2" and not tmpLink in meroXLinksReplicate2:
        meroXLinksReplicate2.append(tmpLink)
    elif row.Replicate == "R3" and not tmpLink in meroXLinksReplicate3:
        meroXLinksReplicate3.append(tmpLink)

print("Total Merox Links= ", str(len(meroX_Links)))
print("Total Merox Links R1= ", str(len(meroXLinksReplicate1)))
print("Total Merox Links R2= ", str(len(meroXLinksReplicate2)))
print("Total Merox Links R3= ", str(len(meroXLinksReplicate3)))


########### Reading all MaxLynx data ################
print("Started: MaxLynx")
maxlynxFile = "C:/yilmaz/Revision/param_scanning/large_scale/pcharm_output/maxlynx_PSMs_0_10_SepFdr_pycharm.txt"
maxlynxFile = "C:/yilmaz/Revision/param_scanning/large_scale/pcharm_output/maxlynx_PSMs_0_10_SepFDR_FDR100_pycharm.txt"


df = pd.read_csv(maxlynxFile, sep=';')
maxlynx_Links, maxlynxLinksReplicate1, maxlynxLinksReplicate2, maxlynxLinksReplicate3 = [], [], [], []

for index, row in df.iterrows():
    tmpReplicate = row.Replicate
    tmpLink = row.LinkedPro1, row.LinkedProSites1, row.LinkedPro2, row.LinkedProSites2
    if not tmpLink in maxlynx_Links:
        maxlynx_Links.append(tmpLink)
    if row.Replicate == "R1" and not tmpLink in maxlynxLinksReplicate1:
        maxlynxLinksReplicate1.append(tmpLink)
    elif row.Replicate == "R2" and not tmpLink in maxlynxLinksReplicate2:
        maxlynxLinksReplicate2.append(tmpLink)
    elif row.Replicate == "R3" and not tmpLink in maxlynxLinksReplicate3:
        maxlynxLinksReplicate3.append(tmpLink)

print("Total MaxLynx Links= ", str(len(maxlynx_Links)))
print("Total MaxLynx Links R1= ", str(len(maxlynxLinksReplicate1)))
print("Total MaxLynx Links R2= ", str(len(maxlynxLinksReplicate2)))
print("Total MaxLynx Links R3= ", str(len(maxlynxLinksReplicate3)))


outputFigure5 = "C:/yilmaz/Revision/Figures/Figure5_MaxLynxVsMeroX_PXD012546.png"
outputFigureSI = "C:/yilmaz/Revision/Figures/SI_MaxLynxVsMeroX_Replicate_PXD012546.png"

outputFigure5 = "C:/yilmaz/Revision/Figures/SI_MaxLynxFDR100_VsMeroX_PXD012546.png"
outputFigureSI = "C:/yilmaz/Revision/Figures/SI_MaxLynxFDR100VsMeroX_Replicate_PXD012546.png"
########################################################################################################
########### Find shared links between MeroX and MaxLynx, when no additional filtering on ###############
print("When no restrictions...")
intersectNoRest = set(meroX_Links).intersection(maxlynx_Links)
print('merox Links = ', str(len(meroX_Links)), "\n", 'maxlynx Links = ', str(len(maxlynx_Links)))
print("intersection Links = " , str(len(list(intersectNoRest))))
onlyMerox = len(meroX_Links) - len(intersectNoRest)
onlyMaxLynx = len(maxlynx_Links) - len(intersectNoRest)
print("only merox Links= " , str(onlyMerox), "\n",  "only maxlynx Links= " , str(onlyMaxLynx), "\n")


allrep = set.union(set.union(set(maxlynxLinksReplicate1), maxlynxLinksReplicate2), maxlynxLinksReplicate3)
print("allrep:" , str(len(allrep)))

fig, ax = plt.subplots(1, 2, figsize=(8,5))

ax[0].text(-0.9, 0.8, "(a)", size=15, weight='bold')
#ax[0].set_title("", pad=-100)
venn3(subsets = (set(maxlynxLinksReplicate1), set(maxlynxLinksReplicate2) , set(maxlynxLinksReplicate3)), ax=ax[0],
      set_labels=('Replicate 1 \n'+ "(" + str(len(maxlynxLinksReplicate1)) + ")", 'Replicate 2 \n'+ "(" + str(len(maxlynxLinksReplicate2)) + ")", 'Replicate 3 \n' + "(" + str(len(maxlynxLinksReplicate3))+ ")"),  alpha = 0.5,
     subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/len(allrep)):1.0%}" + ")")

ax[1].text(-0.6, 0.8, "(b)", size=15, weight='bold')
#ax[1].set_title("MaxLynx vs MeroX")
venn2(subsets = (set(maxlynx_Links), set(meroX_Links)),ax=ax[1],
      set_colors=sns.color_palette("colorblind"),
      set_labels=('MaxLynx \n'+ "(" + str(len(maxlynx_Links)) + ")", 'MeroX \n'+ "(" + str(len(meroX_Links)) + ")"),  alpha = 0.5,
      subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/(onlyMaxLynx+onlyMerox+len(list(intersectNoRest)))):1.0%}" + ")" )

plt.savefig(outputFigure5, dpi=300)
plt.show()

fig, ax = plt.subplots(1, 3, figsize=(10,5))
plt.suptitle("MaxLynx vs MeroX")
ax[0].set_title("no filtering")

venn2(subsets = (set(maxlynx_Links), set(meroX_Links)),ax=ax[0],
      set_colors=sns.color_palette("colorblind"),
      set_labels=('MaxLynx \n'+ "(" + str(len(maxlynx_Links)) + ")", 'MeroX \n'+ "(" + str(len(meroX_Links)) + ")"),  alpha = 0.5,
      subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/(onlyMaxLynx+onlyMerox+len(list(intersectNoRest)))):1.0%}" + ")" )


# ########################################################################################################
# ############ Find shared links between MeroX and MaxLynx, were found in 2/3 replicates #################
print("When reporting cross links existing in 2/3 replicates...")
meroxinreplicate12 = set(meroXLinksReplicate1).intersection(meroXLinksReplicate2)
meroxinreplicate13 = set(meroXLinksReplicate1).intersection(meroXLinksReplicate3)
meroxinreplicate23 = set(meroXLinksReplicate2).intersection(meroXLinksReplicate3)
meroxin23replicates = set(list(meroxinreplicate12) + list(meroxinreplicate13) + list(meroxinreplicate23))

maxlynxinreplicate12 = set(maxlynxLinksReplicate1).intersection(maxlynxLinksReplicate2)
maxlynxinreplicate13 = set(maxlynxLinksReplicate1).intersection(maxlynxLinksReplicate3)
maxlynxinreplicate23 = set(maxlynxLinksReplicate2).intersection(maxlynxLinksReplicate3)
maxlynxin23replicates = set(list(maxlynxinreplicate12) + list(maxlynxinreplicate13) + list(maxlynxinreplicate23))

print('maxlynx in 2/3 replicates= ', str(len(maxlynxin23replicates)))

intersectIn23Replicates = set(meroxin23replicates).intersection(maxlynxin23replicates)
print("intersect-In 2/3 replicates=" , str(len(intersectIn23Replicates)))

onlyMerox = len(meroxin23replicates) - len(intersectIn23Replicates)
onlyMaxLynx = len(maxlynxin23replicates) - len(intersectIn23Replicates)
print("only merox Links= " , str(onlyMerox), "\n", "only maxlynx Links= " , str(onlyMaxLynx))


ax[1].set_title("found in\n 2/3 replicates")
venn2(subsets = (maxlynxin23replicates, meroxin23replicates),ax=ax[1],
      set_colors=sns.color_palette("colorblind"),
      set_labels=('MaxLynx \n'+ "(" + str(len(maxlynxin23replicates)) + ")", 'MeroX \n'+ "(" + str(len(meroxin23replicates)) + ")"),  alpha = 0.5,
      subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/(onlyMaxLynx+onlyMerox+len(intersectIn23Replicates))):1.0%}" + ")" )



########################################################################################################
############ Find shared links between MeroX and MaxLynx, were found in 3/3 replicates ################
print("When reporting cross links existing in 3/3 replicates...")
meroxinallreplicates = set(meroXLinksReplicate3).intersection(meroXLinksReplicate2, meroXLinksReplicate1)
maxlynxinallreplicates = set(maxlynxLinksReplicate1).intersection(maxlynxLinksReplicate2, maxlynxLinksReplicate3)
print('merox in 3/3 replicates = ', str(len(meroxinallreplicates)))
print('maxlynx in 3/3 replicates= ', str(len(maxlynxinallreplicates)))

intersectInAllRep = set(meroxinallreplicates).intersection(maxlynxinallreplicates)
print("intersect-In 3/3 replicates=" , str(len(intersectInAllRep)))

onlyMerox = len(meroxinallreplicates) - len(intersectInAllRep)
onlyMaxLynx = len(maxlynxinallreplicates) - len(intersectInAllRep)
print("only merox Links= " , str(onlyMerox), "\n", "only maxlynx Links= " , str(onlyMaxLynx), "\n")


ax[2].set_title("found in\n 3/3 replicates")
venn2(subsets = (maxlynxinallreplicates, meroxinallreplicates), ax=ax[2],
      set_colors=sns.color_palette("colorblind"),
      set_labels=('MaxLynx \n'+ "(" + str(len(maxlynxinallreplicates)) + ")", 'MeroX \n'+ "(" + str(len(meroxinallreplicates)) + ")"),  alpha = 0.5,
      subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/(onlyMaxLynx+onlyMerox+len(intersectInAllRep))):1.0%}" + ")" )


plt.savefig(outputFigureSI, dpi=300)
plt.show()
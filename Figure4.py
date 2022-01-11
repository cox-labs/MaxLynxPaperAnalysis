import pandas as pd
import matplotlib.pyplot as plt


fdr1crosslinksDsbu4 =  pd.DataFrame({
    "correct":[242, 207, 237, 120], "incorrect":[10, 11, 15, 37]},
    index=["MaxLynx", "MeroX-Rise", "MeroX-Riseup", "XlinkX"])

fdr1crosslinksDsso4 = pd.DataFrame({
    "correct":[185, 124, 149, 128], "incorrect":[3, 1, 19, 53]},
    index=["MaxLynx", "MeroX-Rise", "MeroX-Riseup", "XlinkX"])



fig, ax = plt.subplots(nrows=1, ncols=2,  figsize=(20, 8))

p1 = ax[0].bar( fdr1crosslinksDsbu4.index, fdr1crosslinksDsbu4.iloc[:,0], label=fdr1crosslinksDsbu4.iloc[:,0], color='#2980b9',)
p2 = ax[0].bar( fdr1crosslinksDsbu4.index, fdr1crosslinksDsbu4.iloc[:,1], label=fdr1crosslinksDsbu4.iloc[:,1], bottom= fdr1crosslinksDsbu4.iloc[:,0], color='#e74c3c')

ax[0].set_ylim(0, 275)
ax[0].set_ylabel("#unique cross links", fontsize = 15)
ax[0].set_title("DSBU data set \n ", fontsize = 17)
ax[0].yaxis.grid(False)
ax[0].tick_params(axis='y', which='major', labelsize=13)
ax[0].tick_params(axis='x', which='major', labelsize=15)
ax[0].bar_label(p1, label_type='center', color="black", fontsize=15)
ax[0].bar_label(p2, label_type='center', fontsize=15)
ax[0].text(-0.7, 290 , "(a)", fontsize=20, fontweight = "demibold")


p3 = ax[1].bar( fdr1crosslinksDsso4.index, fdr1crosslinksDsso4.iloc[:,0], label=fdr1crosslinksDsso4.iloc[:,0], color='#2980b9',)
p4 = ax[1].bar( fdr1crosslinksDsso4.index, fdr1crosslinksDsso4.iloc[:,1], label=fdr1crosslinksDsso4.iloc[:,1] , bottom= fdr1crosslinksDsso4.iloc[:,0], color='#e74c3c')

ax[1].set_ylim(0, 275)
ax[1].set_ylabel("#unique cross links", fontsize = 15)
ax[1].legend(["Correct","Incorrect"], fontsize = 15)
ax[1].set_title("DSSO data set \n ", fontsize = 17)
ax[1].yaxis.grid(False)
ax[1].tick_params(axis='y', which='major', labelsize=13)
ax[1].tick_params(axis='x', which='major', labelsize=15)
ax[1].bar_label(p3, label_type='center', color="black", fontsize=15)
ax[1].bar_label(p4, label_type='center', fontsize=15)
ax[1].text(-0.7, 290 , "(b)", fontsize=20, fontweight = "demibold")

plt.savefig('C:/yilmaz/Revision/Figures/Figure4_revised_dsso_dsbu_diff_fdr_001.png', dpi=300)
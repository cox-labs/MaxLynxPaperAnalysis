import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

folder = "C:/yilmaz/Revision/param_scanning/tims-tof/toreport/"
dsbuInput = folder + "tims_dsbu_score_0_partial_10_fdr-separate_crosslinkMsms.txt"
dssoInput = folder + "TIMS_DSSO_score_0_partial_10_fdr-separate_crosslinkMsms.txt"

df = pd.read_csv(dsbuInput, sep='\t')
dfToPlot = df.loc[df["Decoy"] == "forward"]
dfToPlot = dfToPlot.loc[dfToPlot["Crosslink Product Type"] != "InterProXL"]

df2 = pd.read_csv(dssoInput, sep='\t')
dfToPlot2 = df2.loc[df2["Decoy"] == "forward"]
dfToPlot2 = dfToPlot2.loc[dfToPlot2["Crosslink Product Type"] != "InterProXL"]

fig, ax = plt.subplots(2, 1, figsize=(13, 20), sharex=True)
sns.set_context("talk", font_scale=1)

ax[0].text(-500, 1500, "(a)", size=23, weight='bold')
ax[0] = sns.scatterplot(data=dfToPlot, x="Mass", y="CCS", hue="Crosslink Product Type", style="Charge", s=120,
                        ax=ax[0]);
ax[0].set_xlim(0, 6000)
ax[0].set_ylim(0, 1400)
ax[0].set_title("DSBU K-KSTY")
ax[0].set_ylabel("CCS (Å2)")
ax[0].legend(bbox_to_anchor=(1.35, 1), borderaxespad=0)

ax[1].text(-500, 1500, "(b)", size=23, weight='bold')
ax[1] = sns.scatterplot(data=dfToPlot, x="Mass", y="CCS", hue="Crosslink Product Type", style="Charge", s=120,
                        ax=ax[1]);
ax[1].set_xlim(0, 6000)
ax[1].set_ylim(0, 1400)
ax[1].set_title("DSSO K-KSTY")
ax[1].set_xlabel("Mass (Da)")
ax[1].set_ylabel("CCS (Å2)")
legend = ax[1].legend(bbox_to_anchor=(1.35, 1), borderaxespad=0)
legend.remove()

fig.savefig("C:/yilmaz/Revision/Figures/Figure6_Timstof.png", format='png', dpi=300, bbox_inches='tight')
plt.show()

import os
import pandas as pd


fastGear_out = pd.read_csv(snakemake.input[0],sep="\s+",header=1)
fastGear_out = fastGear_out[["StrainName", "Start", "End"]]

fastGear_out.to_csv(snakemake.output[0], header=False, index=False, sep="\t")





    
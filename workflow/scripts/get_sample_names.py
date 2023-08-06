import os

with open(snakemake.output[0], "w") as f:
    for file in os.listdir(snakemake.input[0]):
        f.write(os.path.splitext(file)[0] + "\n")

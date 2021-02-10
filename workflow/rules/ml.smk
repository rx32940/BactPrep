rule iqtree:
    input:
        msa = rules.roary.output if config["within_species"] else rules.pirate.output
    output:
        iqtree_dir + project_prefix +".treefile"
    conda:
        "../env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix = iqtree_dir + project_prefix
    shell:
        """
        iqtree -nt AUTO -m MFP -pre {params.prefix} -s {input}
        """

    
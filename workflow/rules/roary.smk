rule Roary:
    input:
        expand(rules.get_gff.output, sample=SAMPLES)
    output:
        roary_dir + "core_gene_alignment.aln",
        roary_dir + "core_alignment_header.embl",
        roary_dir + "gene_presence_absence.csv"
    threads:
        THREADS
    conda:
        WORKFLOW + "env/roary.yaml"
    params:
        alignment = "-e -n",
        out_dir = roary_dir,
        additional=" " + roary_params if roary_params != "" else ""
    shell:
        """
            roary{params.additional} -e -n -p {threads} -f {params.out_dir} {input} -v -cd {core_percentage} -z

            cp -r {roary_dir}*/* {roary_dir}

            rm -rf {roary_dir}_*
        """

rule core_gene_concatenation_ML_tree:
    input:
        msa = rules.Roary.output[0]
    output:
        iqtree_dir + project_prefix +".treefile"
    conda:
        WORKFLOW + "env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix = iqtree_dir + project_prefix
    shell:
        """
        iqtree -bb 1000 -nt AUTO -m MFP -pre {params.prefix} -s {input}
        """

    


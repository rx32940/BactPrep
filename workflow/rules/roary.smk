rule Roary:
    input:
        expand(rules.get_gff.output, sample=SAMPLES)
    output:
        roary_dir + "core_gene_alignment.aln"
    threads:
        THREADS
    conda:
        "../env/roary.yaml"
    params:
        alignment = "-e -n",
        out_dir = roary_dir
    shell:
        """
            roary -e -n -p {threads} -f {params.out_dir} {input} -v -cd {core_percentage} -z

            cp -r {roary_dir}*/* {roary_dir}
        """

rule core_gene_concatenation_ML_tree:
    input:
        msa = rules.Roary.output
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
        iqtree -bb 1000 -nt AUTO -m MFP -pre {params.prefix} -s {input}
        """

    


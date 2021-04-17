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

    


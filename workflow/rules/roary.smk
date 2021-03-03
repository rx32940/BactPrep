rule roary:
    input:
        rules.get_gff.output
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
            roary -e -n -p {threads} -f {params.out_dir} {input}/*.gff -v -cd {core_percentage}

            cp -r {roary_dir}*/* {roary_dir}
        """

    


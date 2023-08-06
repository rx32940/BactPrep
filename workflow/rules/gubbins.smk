rule Gubbins:
    input:
        rules.clean_snippy_core.output
    output:
        gubbins_dir + project_prefix + ".filtered_polymorphic_sites.fasta"
    conda:
        WORKFLOW + "env/gubbins.yaml"
    threads:
        THREADS
    params:
        out_prefix = gubbins_dir + project_prefix,
        additional=" " + gubbins_params if gubbins_params != "" else ""
    shell:
        """
        cd {gubbins_dir}
        run_gubbins.py{params.additional} --threads {threads} -v -p {params.out_prefix} {input}
        """

rule clean_snps:
    input:
        gubbins_aln=rules.Gubbins.output
    output:
        gubbins_dir + ".filtered_polymorphic_sites.noRef.fasta"
    conda:
        WORKFLOW + "env/gubbins.yaml"
    shell:
        """
        cat {input} | seqkit grep -v -p Reference > {output}
        """

rule Gubbins_SNPS_ML_tree:
    input:
        rules.clean_snps.output
    output:
        gubbins_dir + "iqtree/" + project_prefix + ".gubbins.treefile"
    conda:
        WORKFLOW + "env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix=gubbins_dir + "iqtree/" + project_prefix + ".gubbins"
    shell:
        """
        iqtree -bb 1000 -nt AUTO -m MFP -pre {params.prefix} -s {input}
        """
    


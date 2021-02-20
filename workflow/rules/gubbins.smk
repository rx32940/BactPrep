rule gubbins:
    input:
        rules.clean_snippy_core.output
    output:
        gubbins_dir + project_prefix + ".filtered_polymorphic_sites.fasta"
    conda:
        "../env/gubbins.yaml"
    threads:
        THREADS
    params:
        out_prefix = gubbins_dir + project_prefix
    shell:
        """
        run_gubbins.py --threads {threads} -v -p {params.out_prefix} {input}
        """

rule clean_snps:
    input:
        rules.gubbins.output
    output:
        gubbins_dir + "snp-sites/" + project_prefix + ".recombFreeSnpsAtcg.fasta"
    conda:
        "../env/gubbins.yaml"
    shell:
        """
        cat {input} | seqkit grep -v -p Reference > {gubbins_dir}{project_prefix}_noref.filtered_polymorphic_sites.fasta
        snp-sites -c {gubbins_dir}{project_prefix}_noref.filtered_polymorphic_sites.fasta > {output}
        """

rule SNPS_tree:
    input:
        rules.clean_snps.output
    output:
        gubbins_dir + "iqtree/" + project_prefix + ".recombFreeSnpsAtcg.treefile"
    conda:
        "../env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix=gubbins_dir + "iqtree/" + project_prefix + ".recombFreeSnpsAtcg"
    shell:
        """
        iqtree -nt AUTO -m MFP+ASC -pre {params.prefix} -s {input}
        """
    


        

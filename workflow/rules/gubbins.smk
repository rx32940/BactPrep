rule Gubbins:
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
        cd {gubbins_dir}
        run_gubbins.py --threads {threads} -v -p {params.out_prefix} {input}
        """

rule clean_snps:
    input:
        rules.Gubbins.output
    output:
        gubbins_dir + "snp-sites/" + project_prefix + ".recombFreeSnpsAtcg.fasta"
    conda:
        "../env/gubbins.yaml"
    shell:
        """
        cat {input} | seqkit grep -v -p Reference > {gubbins_dir}{project_prefix}_noref.filtered_polymorphic_sites.fasta
        snp-sites -c {gubbins_dir}{project_prefix}_noref.filtered_polymorphic_sites.fasta > {output}
        """

rule Gubbins_SNPS_ML_tree:
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
    

rule annotate_GubbinsSNP_alignment:
    input:
        original_alignment = rules.clean_snps.output, 
        metadata_file = metadata_file
    output:
        gubbins_dir  + project_prefix + "_meta.recombFreeSnpsAtcg.fasta"
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1
    shell:
        """
        python {workflow.basedir}/scripts/change_fasta_header.py \
        {input.metadata_file} {input.original_alignment} {params.meta_include} {params.key_column_index} {output}
        """

rule annotate_GubbinsSNPs_tree:
    input:
        tree=rules.Gubbins_SNPS_ML_tree.output,
        metadata_file = metadata_file
    output:
        gubbins_dir  + "iqtree/" +  project_prefix + "_meta.GubbinsSNPs.newick"
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1
    shell:
        """
        python {workflow.basedir}/scripts/rename_phylogeny_taxa.py \
        {input.metadata_file} {input.tree} {params.meta_include} {params.key_column_index} {output}
        """








        

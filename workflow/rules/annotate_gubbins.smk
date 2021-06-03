rule annotate_GubbinsSNP_alignment:
    input:
        original_alignment = rules.clean_snps.output, 
        metadata_file = metadata_file
    output:
        gubbins_dir  + project_prefix + "_meta.recombFreeSnpsAtcg.fasta"
    params:
        meta_include = metadata_include if metadata_include is not None else None,
        key_column_index = biosample_column -1 if biosample_column is not None else None
    shell:
        """
        if metadata_file != "":
            python {workflow.basedir}/scripts/change_fasta_header.py \
            {input.metadata_file} {input.original_alignment} {params.meta_include} {params.key_column_index} {output}
        else:
            touch {output}
        """

rule annotate_GubbinsSNPs_tree:
    input:
        tree=rules.Gubbins_SNPS_ML_tree.output,
        metadata_file = metadata_file
    output:
        gubbins_dir  + "iqtree/" +  project_prefix + "_meta.GubbinsSNPs.newick"
    params:
        meta_include = metadata_include if metadata_include is not None else None,
        key_column_index = biosample_column -1 if biosample_column is not None else None
    shell:
        """
        if metadata_file != "":
            python {workflow.basedir}/scripts/rename_phylogeny_taxa.py \
            {input.metadata_file} {input.tree} {params.meta_include} {params.key_column_index} {output}
        else:
            touch {output}
        """









        

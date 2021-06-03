rule annotate_coreGene_tree:
    input:
        tree=rules.core_gene_concatenation_ML_tree.output,
        metadata_file = metadata_file
    output:
        iqtree_dir + project_prefix + "_meta.coreConcate.newick"
    params:
        meta_include = metadata_include if metadata_include is not None else None,
        key_column_index = biosample_column -1 if metadata_include is not None else None
    shell:
        """
        if metadata_file != "":
            python {workflow.basedir}/scripts/rename_phylogeny_taxa.py \
            {input.metadata_file} {input.tree} {params.meta_include} {params.key_column_index} {output}
        else:
            touch {output}
        """
rule core_gene_concatenation_ML_tree:
    input:
        msa = rules.Roary.output if config["within_species"] else rules.PIRATE.output
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

rule annotate_coreGene_tree:
    input:
        tree=rules.core_gene_concatenation_ML_tree.output,
        metadata_file = metadata_file
    output:
        iqtree_dir + project_prefix + "_meta.coreSNPs.newick"
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1
    shell:
        """
        python {workflow.basedir}/scripts/rename_phylogeny_taxa.py \
        {input.metadata_file} {input.tree} {params.meta_include} {params.key_column_index} {output}
        """

    
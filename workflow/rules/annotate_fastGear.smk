rule annotate_core_SNPs_alignment:
    input:
        original_alignment = rules.call_snp_from_masked_alignment.output,
        metadata_file = metadata_file
    output:
        os.path.join(fastGear_dir , str(project_prefix + "_core_mask_snp_meta.fasta"))
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1
    shell:
        """
        if [[ -z {input.metadata_file} ]]; 
        then
            python {workflow.basedir}/scripts/change_fasta_header.py {input.metadata_file} {input.original_alignment} {params.meta_include} {params.key_column_index} {output}
        else
            touch {output}
        fi
        """

rule annotate_coreSNP_tree:
    input:
        tree=rules.core_genome_snps_ML_tree.output,
        metadata_file = metadata_file
    output:
        os.path.join(fastGear_dir , "fastgear_iqtree" , str(project_prefix + "_meta.coreSNPs.newick"))
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1
    shell:
        """
        if [[ -n {input.metadata_file} ]]; 
        then
            python {workflow.basedir}/scripts/rename_phylogeny_taxa.py {input.metadata_file} {input.tree} {params.meta_include} {params.key_column_index} {output}
        else
            touch {output}
        fi
        """
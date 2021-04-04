rule fastGear:
    input:
        rules.Roary.output
    output:
        fastGear_dir + "output/recombinations_recent.txt"
    params:
        output_mat = fastGear_dir + project_prefix + ".mat"
    shell:
        """
        LD_LIBRARY_PATH={matlab_path}
        {fastGear_exe}run_fastGEAR.sh {mcr_path} {input} {params.output_mat} {fastGear_exe}fG_input_specs.txt 
        """

rule convert_recombination_to_bed:
    input:
        recent = rules.fastGear.output
    output:
        fastGear_dir + "output/recombinations_recent.bed"
    script:
        "../scripts/fastGear_to_bed.py"
 

rule mask_recombination:
    input:
        bed = rules.convert_recombination_to_bed.output,
        fasta = rules.Roary.output
    output:
        fastGear_dir +  project_prefix + "_core_mask.fasta"
    conda:
        "../env/bedtools.yaml"
    shell:
        """
        bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}
        """

rule call_snp_from_masked_alignment:
    input:
        rules.mask_recombination.output
    output:
        fastGear_dir +  project_prefix + "_core_mask_snp.fasta"
    conda:
        "../env/bedtools.yaml"
    shell:
        """
            snp-sites -c {input} -o {output}
        """

rule core_genome_snps_ML_tree:
    input:
        rules.call_snp_from_masked_alignment.output
    output:
        fastGear_dir + "fastgear_iqtree/" + project_prefix + "_core_mask_snp.treefile"
    conda:
        "../env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix = fastGear_dir + "fastgear_iqtree/" + project_prefix + "_core_mask_snp"
    shell:
        """
            iqtree -nt AUTO -m MFP+ASC -pre {params.prefix} -s {input}
        """

rule annotate_core_SNPs_alignment:
    input:
        original_alignment = rules.call_snp_from_masked_alignment.output,
        metadata_file = metadata_file
    output:
        fastGear_dir + project_prefix + "_core_mask_snp_meta.fasta"
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1  
    shell:
        """
            python {workflow.basedir}/scripts/change_fasta_header.py \
            {input.metadata_file} {input.original_alignment} {params.meta_include} {params.key_column_index} {output}
        """

rule annotate_coreSNP_tree:
    input:
        tree=rules.core_genome_snps_ML_tree.output,
        metadata_file = metadata_file
    output:
        fastGear_dir + "fastgear_iqtree/" +  project_prefix + "_meta.coreSNPs.newick"
    params:
        meta_include = metadata_include,
        key_column_index = biosample_column -1
    shell:
        """
        python {workflow.basedir}/scripts/rename_phylogeny_taxa.py \
        {input.metadata_file} {input.tree} {params.meta_include} {params.key_column_index} {output}
        """






        

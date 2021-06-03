rule fastGear:
    input:
        rules.Roary.output
    output:
        os.path.join(fastGear_dir , "output/recombinations_recent.txt")
    params:
        output_mat = os.path.join(fastGear_dir , str(project_prefix + ".mat"))
    shell:
        """
        LD_LIBRARY_PATH={matlab_path}
        {fastGear_exe}run_fastGEAR.sh {mcr_path} {input} {params.output_mat} {fastGear_exe}fG_input_specs.txt 
        """

rule convert_recombination_to_bed:
    input:
        recent = rules.fastGear.output
    output:
        os.path.join(fastGear_dir , "output/recombinations_recent.bed")
    script:
        "../scripts/fastGear_to_bed.py"
 

rule mask_recombination:
    input:
        bed = rules.convert_recombination_to_bed.output,
        fasta = rules.Roary.output
    output:
        os.path.join(fastGear_dir , str(project_prefix + "_core_mask.fasta"))
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
        os.path.join(fastGear_dir , str(project_prefix + "_core_mask_snp.fasta"))
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
        os.path.join(fastGear_dir , "fastgear_iqtree" , str(project_prefix + "_core_mask_snp.treefile"))
    conda:
        "../env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix = os.path.join(fastGear_dir , "fastgear_iqtree" , str(project_prefix + "_core_mask_snp"))
    shell:
        """
            iqtree -bb 1000 -nt AUTO -m MFP+ASC -pre {params.prefix} -s {input}
        """








        

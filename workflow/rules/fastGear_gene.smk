rule move_input_directory:
    output: 
        os.path.join(fastGear_gene_dir, "input_alns", "{gene}.fasta")
    params:
        input_files=config["fastgear_gene_file_list"]
    shell:
        """
        for i in {params.input_files};
        do
        cp $i {fastGear_gene_dir}input_alns/{wildcards.gene}.fasta
        done
        """

rule fastGear_singleGene:
    input:
        rules.move_input_directory.output
    output:
        os.path.join(fastGear_gene_dir, "{gene}", str("{gene}.mat"))
    params:
        output_mat = os.path.join(fastGear_gene_dir, "{wildcards.gene}", str("{wildcards.gene}.mat")),
        input_files = config["fastgear_gene_file_list"]
    shell:
        """
        LD_LIBRARY_PATH={matlab_path}
        {fastGear_exe}run_fastGEAR.sh {mcr_path} {input} {fastGear_gene_dir}{wildcards.gene}/{wildcards.gene}.mat {fastGear_exe}fG_input_specs.txt 
        """

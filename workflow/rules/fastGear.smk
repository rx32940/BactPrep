rule fastGear:
    input:
        rules.roary.output
    output:
        fastGear_dir + "roary_recomb_detect.mat"
    shell:
        """
        
        LD_LIBRARY_PATH={matlab_path}

        {fastGear_exe}run_fastGEAR.sh {mcr_path} {input} {output} {fastGear_exe}fG_input_specs.txt 

        """
rule PIRATE:
    input:
        gff_dir
    output: 
        os.path.join(pirate_dir , "core_alignment.fasta")
    conda: 
        "../env/pirate.yaml"
    threads: 
        THREADS
    params:
        out_dir = pirate_dir,
        mcl="-pan-opt \"-f 6\""
    shell:
        """
            PIRATE \
            -i {input} -o {params.out_dir} \
            -a -r -t {threads} \
            {params.mcl}       
        """

rule convert_to_Roary:
    input:
        os.path.join(pirate_dir , "PIRATE.*.tsv")
    output:
        os.path.join(pirate_dir , "pirate_roary_pres_abs.csv")
    params:
        script_dir=workflow.basedir + "/scripts/"
    shell:
        """
           perl {params.script_dir}PIRATE_to_roary.pl -i {input} -o {pirate_dir}pirate_roary_pres_abs.csv 
        
        """
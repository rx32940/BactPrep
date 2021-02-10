rule pirate:
    input:
        gff_dir
    output: 
        pirate_dir + "core_alignment.fasta"
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

rule convert_roary:
    input:
        pirate_dir + "PIRATE.*.tsv"
    output:
        pirate_dir + "pirate_roary_pres_abs.csv"
    params:
        script_dir=workflow.basedir + "/scripts/"
    shell:
        """
           perl {params.script_dir}PIRATE_to_roary.pl -i {input} -o {pirate_dir}pirate_roary_pres_abs.csv 
        
        """
rule pirate:
    input: 
        expand([prokka_dir + "{sample}/{sample}.gff"], sample= SAMPLES)
    output: 
        pirate_dir + "pirate_roary_pres_abs.csv"
    conda: 
        "../env/pirate.yaml"
    threads: 
        THREADS
    params:
        out_dir = pirate_dir
        mcl="-pan-opt \"-f 6\""
        script_dir="../scripts/"
    shell:
        """
            PIRATE \
            -i {input} -o {params.out_dir} \
            -a -r -t {threads} \
            {params.mcl} 

            perl {script_dir}PIRATE_to_roary.pl -i {pirate_dir}PIRATE.*.tsv -o {pirate_dir}pirate_roary_pres_abs.csv
        """
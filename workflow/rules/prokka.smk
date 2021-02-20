if config["prokka_dir"]:
    rule get_gff:
        output:
            directory(gff_dir)
        shell:
            """
            mkdir {gff_dir}
            cp {prokka_dir}*/*.gff {output}
            """

else:
    rule prokka:
        input:
            asm_dir + "/{sample}.fna"
        output:
            directory(prokka_dir + "{sample}")
        conda:
            "../env/prokka.yaml"
        threads:
            THREADS
        params:
            genus = config["genus"],
            kingdom = config["kingdom"],
            output_dir = prokka_dir
        shell:
            """
            prokka -kingdom {params.kingdom} -genus {params.genus} \
            -outdir {output} \
            -prefix {wildcards.sample} \
            {input} --force -cpu {threads}
            
            """
    rule get_gff:
        output:
            directory(gff_dir)
        shell:
            """
                mkdir {gff_dir}
                cp {prokka_dir}*/*.gff {output}
            """


        

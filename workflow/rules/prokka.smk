rule rename_file_prefix:
    output:
        os.path.join(asm_dir , "{sample}.fna")
    shell:
        """
        mv {asm_dir}{wildcards.sample}* {output}
        """

rule Prokka:
    input:
        os.path.join(asm_dir , "{sample}.fna")
    output:
        directory(os.path.join(prokka_dir , "{sample}"))
    conda:
        WORKFLOW + "env/prokka.yaml"
    threads:
        THREADS
    params:
        kingdom = config["kingdom"],
        output_dir = prokka_dir
    shell:
        """
        prokka -kingdom {params.kingdom} \
        -outdir {output} \
        -prefix {wildcards.sample} \
        {input} --force -cpu {threads}
        
        """
rule get_gff:
    input:
        rules.Prokka.output
    output:
        os.path.join(gff_dir , "{sample}.gff")
    shell:
        """
            cp {input}/{wildcards.sample}.gff {gff_dir}{wildcards.sample}".gff"
        """


    

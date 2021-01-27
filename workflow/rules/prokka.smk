rule prokka:
    input:
        asm_dir + "/{sample}.fna"
    output:
        prokka_dir + "{sample}/{sample}.gff"
    conda:
        "../env/prokka.yaml"
    params:
        genus = config["genus"],
        kingdom = config["kingdom"],
        output_dir = prokka_dir
    shell:
        """
        prokka -kingdom {params.kingdom} -genus {params.genus} \
        -outdir {params.output_dir}{wildcards.sample} \
        -prefix {wildcards.sample} \
        {input} --force
        """


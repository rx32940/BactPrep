
rule snippy_multi:
    input:
        os.path.join(asm_dir , "{sample}.fna")
    output:
        directory(os.path.join(snippy_dir , "{sample}"))
    conda:
        "../env/snippy.yaml"
    threads:
        THREADS
    shell:
        """
        snippy --outdir {output} --ctgs {input} --ref {reference} --cpus {threads}
        """

rule snippy_core:
    input:
        expand([snippy_dir + "{sample}"], sample = SAMPLES)
    output:
        os.path.join(snippy_dir , "core.full.aln")
    conda:
        "../env/snippy.yaml"
    params:
        mask=str("--mask ") + phage if phage is not None else None
    shell:
        """
        snippy-core {params.mask} --ref {reference} {input} --prefix {snippy_dir}core
        """

rule clean_snippy_core:
    input:
        os.path.join(snippy_dir , "core.full.aln")
    output:
        os.path.join(snippy_dir , "clean.full.aln")
    conda:
        "../env/snippy.yaml"
    shell:
        """
        snippy-clean_full_aln {input} > {snippy_dir}clean.full.aln
        
        """








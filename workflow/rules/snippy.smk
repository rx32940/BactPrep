
rule snippy_multi:
    input:
        asm_dir + "/{sample}.fna"
    output:
        directory(snippy_dir + "{sample}")
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
        snippy_dir + "core.full.aln"
    conda:
        "../env/snippy.yaml"
    shell:
        """
        snippy-core --mask {phage} --ref {reference} {input} --prefix {snippy_dir}core
        """

rule clean_snippy_core:
    input:
        snippy_dir + "core.full.aln"
    output:
        snippy_dir + "clean.full.noref.aln"
    conda:
        "../env/snippy.yaml"
    shell:
        """
        snippy-clean_full_aln {input} > {snippy_dir}clean.full.aln
        
        """








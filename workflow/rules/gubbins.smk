rule Gubbins:
    input:
        rules.clean_snippy_core.output
    output:
        gubbins_dir + project_prefix + ".recombination_predictions.gff"
    conda:
        "../env/gubbins.yaml"
    threads:
        THREADS
    params:
        out_prefix = gubbins_dir + project_prefix,
        additional=" " + gubbins_params if gubbins_params != "" else ""
    shell:
        """
        cd {gubbins_dir}
        run_gubbins.py{params.additional} --threads {threads} -v -p {params.out_prefix} {input}
        """

rule clean_snps:
    input:
        gubbins_gff=rules.Gubbins.output,
        snippy_aln=rules.clean_snippy_core.output
    output:
        snippy_dir + "clean.full.recomMaked.noRef.aln"
    conda:
        "../env/gubbins.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/mask_gubbins_aln.py --aln {input.snippy_aln} \
        --gff {input.gubbins_gff} \
        --out {snippy_dir}clean.full.recomMaked.aln

        cat {snippy_dir}clean.full.recomMaked.aln | seqkit grep -v -p Reference > {output}
        """

rule Gubbins_SNPS_ML_tree:
    input:
        rules.clean_snps.output
    output:
        gubbins_dir + "iqtree/" + project_prefix + ".recombFreeSnpsAtcg.treefile"
    conda:
        "../env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix=gubbins_dir + "iqtree/" + project_prefix + ".recombFreeSnpsAtcg"
    shell:
        """
        snp-sites -c {input} -o {gubbins_dir}clean.recomMaked.noRef.SNPs.atcg.aln
        
        iqtree -bb 1000 -nt AUTO -m MFP -pre {params.prefix} -s {gubbins_dir}clean.recomMaked.noRef.SNPs.atcg.aln -fconst $(snp-sites -C {input})
        """
    


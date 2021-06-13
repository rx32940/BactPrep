checkpoint extract_core_loci:
    input:
        rules.Roary.output[1]
    output:
        os.path.join(fastGear_core_dir, "roary_coreGeneAln_locustag.txt")
    shell:
        """
        python {workflow.basedir}/scripts/get.genes.roary.core.aln.py {input} {output}
        """


rule fastGear_core:
    input:
        rules.change_pan_gene_aln_headers.output
    output:
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/output/recombinations_recent.txt"),
        temp(os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/{core_locus}.fa"))
    shell:
        """
        LD_LIBRARY_PATH={matlab_path}
        {fastGear_exe}run_fastGEAR.sh {mcr_path} {fastGear_dir}roary_pangenome_seq/{wildcards.core_locus}.fa.aln {fastGear_core_dir}core_loci_fastGear_out/{wildcards.core_locus}/{wildcards.core_locus}.mat {fastGear_params}
        cp {fastGear_dir}roary_pangenome_seq/{wildcards.core_locus}.fa.aln {fastGear_core_dir}core_loci_fastGear_out/{wildcards.core_locus}/{wildcards.core_locus}.fa
        """


rule plot_core_fastGear:
    input:
        loci=rules.extract_core_loci.output,
        tree=rules.core_gene_concatenation_ML_tree.output,
        fastGear_out=rules.fastGear_core.output
    output:
        fastGear_core_dir + "plot_coregenome/core_fastgear_plot_recombination_count.pdf"
    shell:
        """
        cd {fastGear_core_dir}plot_coregenome/
        python {workflow.basedir}/scripts/post_fastGear.py \
        -i {fastGear_core_dir}core_loci_fastGear_out \
        -g {input.loci} \
        -o {fastGear_core_dir}plot_coregenome/core_fastgear_plot \
        -s True -f pdf -p {input.tree} -z True -y 100 -x 100
        """

rule convert_recombination_to_bed:
    input:
        recent = rules.fastGear_core.output[0]
    output:
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/output/recombinations_recent.bed")
    script:
        "../scripts/fastGear_to_bed.py"
        
rule mask_recombination:
    input:
        bed = rules.convert_recombination_to_bed.output,
        fasta = rules.fastGear_core.output[1]
    output:
        temp(os.path.join(fastGear_core_dir , "masked_coregene_aln","{core_locus}_core_mask.fasta"))
    conda:
        "../env/bedtools.yaml"
    shell:
        """
        bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}
        """

def get_masked_core_aln(wildcards):
    with open(checkpoints.extract_core_loci.get().output[0]) as f:
        core_loci=[locus for locus in f.read().split('\n') if len(locus) > 0]
    return expand(os.path.join(fastGear_core_dir , "masked_coregene_aln","{core_locus}_core_mask.fasta"), core_locus=core_loci)

rule concate_gene_aln:
    input:
        get_masked_core_aln
    output:
        os.path.join(fastGear_core_dir, "fastGear_masked_coregene_aln.fasta")
    shell:
        """
            touch {output}
            cat {input[0]}
            while read $file;
            do
            cat $file >> {output}
            done
        """

rule call_snp_from_masked_alignment:
    input:
        rules.concate_gene_aln.output
    output:
        os.path.join(fastGear_core_dir , str(project_prefix + "_core_mask_snp.fasta"))
    conda:
        "../env/bedtools.yaml"
    shell:
        """
            snp-sites -c {input} -o {output}
        """

rule core_genome_snps_ML_tree:
    input:
        rules.call_snp_from_masked_alignment.output
    output:
        os.path.join(fastGear_core_dir , "fastgear_iqtree" , str(project_prefix + "_core_mask_snp.treefile"))
    conda:
        "../env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix = os.path.join(fastGear_core_dir , "fastgear_iqtree" , str(project_prefix + "_core_mask_snp"))
    shell:
        """
            iqtree -bb 1000 -nt AUTO -m MFP+ASC -pre {params.prefix} -s {input}
        """



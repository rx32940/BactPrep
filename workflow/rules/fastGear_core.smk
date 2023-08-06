checkpoint extract_core_loci:
    input:
        rules.Roary.output[1]
    output:
        os.path.join(fastGear_core_dir, "roary_coreGeneAln_locustag.txt")
    shell:
        """
        python {WORKFLOW}scripts/get.genes.roary.core.aln.py {input} {output}
        """


rule fastGear_core:
    input:
        rules.change_pan_gene_aln_headers.output
    output:
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/output/recombinations_recent.txt"),
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/output/recombinations_ancestral.txt"),
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/output/lineage_information.txt"),
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/{core_locus}.fa")
    shell:
        """
        LD_LIBRARY_PATH={matlab_path}
        {fastGear_exe}run_fastGEAR.sh {mcr_path} {fastGear_dir}roary_pangenome_seq/{wildcards.core_locus}.fa.aln {fastGear_core_dir}core_loci_fastGear_out/{wildcards.core_locus}/{wildcards.core_locus}.mat {fastGear_params}
        cp {fastGear_dir}roary_pangenome_seq/{wildcards.core_locus}.fa.aln {fastGear_core_dir}core_loci_fastGear_out/{wildcards.core_locus}/{wildcards.core_locus}.fa
        """

rule convert_recombination_to_bed:
    input:
        recent = rules.fastGear_core.output[0],
        ancestral=rules.fastGear_core.output[1],
        lineage=rules.fastGear_core.output[2]
    output:
        os.path.join(fastGear_core_dir , "core_loci_fastGear_out","{core_locus}/output/recombinations.bed")
    script:
        WORKFLOW + "scripts/fastGear_to_bed.py"
        
rule mask_recombination:
    input:
        bed = rules.convert_recombination_to_bed.output,
        fasta = rules.fastGear_core.output[3]
    output:
        os.path.join(fastGear_core_dir , "masked_coregene_aln","{core_locus}_core_mask.fasta")
    conda:
        WORKFLOW + "env/bedtools.yaml"
    shell:
        """
        bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}
        """

def get_masked_core_aln(wildcards):
    with open(checkpoints.extract_core_loci.get().output[0]) as f:
        core_loci=[locus for locus in f.read().split('\n') if len(locus) > 0]
    return expand(os.path.join(fastGear_core_dir , "masked_coregene_aln","{core_locus}_core_mask.fasta"), core_locus=core_loci)


rule create_core_loci_file:
    input:
        get_masked_core_aln
    output:
        fastGear_core_dir + "core_loci_list.txt"
    run:
        loci = [file.split("/")[-1].split("_core_mask")[0].strip() for file in input if not file.startswith(".")]
        with open(str(output[0]), mode='w', encoding='utf-8') as myfile:
            myfile.write('\n'.join(loci))


rule plot_core_fastGear:
    input:
        loci=rules.create_core_loci_file.output,
        tree=rules.core_gene_concatenation_ML_tree.output,
    output:
        fastGear_core_dir + "plot_coregenome/core_fastgear_plot_recombination_count.pdf"
    shell:
        """
        cd {fastGear_core_dir}plot_coregenome/
        python {WORKFLOW}scripts/post_fastGear.py \
        -i {fastGear_core_dir}core_loci_fastGear_out \
        -g {input.loci} \
        -o {fastGear_core_dir}plot_coregenome/core_fastgear_plot \
        -s True -f pdf -p {input.tree} -z True -y 100 -x 100
        """

rule concate_gene_aln:
    input:
        get_masked_core_aln
    output:
        os.path.join(fastGear_core_dir, "fastGear_masked_coregene_aln.fasta")
    shell:
        """
            perl {WORKFLOW}scripts/catfasta2phyml.pl -f -c -s --verbose {fastGear_core_dir}masked_coregene_aln/*fasta > {output}
        """

rule call_snp_from_masked_alignment:
    input:
        rules.concate_gene_aln.output
    output:
        os.path.join(fastGear_core_dir , str(project_prefix + "_core_mask_snp.fasta"))
    conda:
        WORKFLOW + "env/bedtools.yaml"
    shell:
        """
            snp-sites {input} -o {output}
        """

rule core_genome_snps_ML_tree:
    input:
        snps_aln=rules.call_snp_from_masked_alignment.output,
        core_aln=rules.concate_gene_aln.output
    output:
        os.path.join(fastGear_core_dir , "fastgear_iqtree" , str(project_prefix + "_core_mask_snp.treefile"))
    conda:
        WORKFLOW + "env/iqtree.yaml"
    threads:
        THREADS
    params:
        prefix = os.path.join(fastGear_core_dir , "fastgear_iqtree" , str(project_prefix + "_core_mask_snp"))
    shell:
        """
            iqtree -bb 1000 -nt AUTO -m MFP -pre {params.prefix} -s {input.snps_aln} -fconst $(snp-sites -C {input.core_aln})
        """



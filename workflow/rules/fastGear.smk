checkpoint change_pan_gene_aln_headers:
    input:
        rules.Roary.output[2]
    output:
        directory(os.path.join(fastGear_dir ,"roary_pangenome_seq"))
    shell:
        """
            mkdir -p {fastGear_dir}roary_pangenome_seq
            python {workflow.basedir}/scripts/change.roary.gene.alns.headers.py {roary_dir}pan_genome_sequences {fastGear_dir}roary_pangenome_seq {input}
        """

rule fastGear:
    input:
        rules.change_pan_gene_aln_headers.output
    output:
        os.path.join(fastGear_dir , "loci_fastGear_out","{locus}/output/recombinations_recent.txt"),
        os.path.join(fastGear_dir , "loci_fastGear_out","{locus}/{locus}.fa")
    shell:
        """
        LD_LIBRARY_PATH={matlab_path}
        {fastGear_exe}run_fastGEAR.sh {mcr_path} {fastGear_dir}roary_pangenome_seq/{wildcards.locus}.fa.aln {fastGear_dir}loci_fastGear_out/{wildcards.locus}/{wildcards.locus}.mat {fastGear_params} 
        cp {fastGear_dir}roary_pangenome_seq/{wildcards.locus}.fa.aln {fastGear_dir}loci_fastGear_out/{wildcards.locus}/{wildcards.locus}.fa
        """

def get_all_loci(wildcards):
    all_loci = [locus.replace("/",".").split(".")[-3] for locus in os.listdir(checkpoints.change_pan_gene_aln_headers.get().output[0]) if not locus.startswith('.') ]
    print(all_loci)
    return expand(os.path.join(fastGear_dir , "loci_fastGear_out","{locus}/output/recombinations_recent.txt"), locus=all_loci)

rule create_loci_file:
    input:
        get_all_loci
    output:
        fastGear_dir + "plot_pangenome/all_loci_list.txt"
    run:
        loci = [file.split("/")[-3] for file in input]
        print(loci)
        with open(str(output[0]), mode='w', encoding='utf-8') as myfile:
            myfile.write('\n'.join(loci))

rule plot_pan_fastGear:
    input:
        loci=rules.create_loci_file.output,
        tree=rules.core_gene_concatenation_ML_tree.output
    output:
        fastGear_dir + "plot_pangenome/pan_fastgear_plot_recombination_count.pdf"

    shell:
        """
        cd {fastGear_dir}plot_pangenome/
        python {workflow.basedir}/scripts/post_fastGear.py \
        -i {fastGear_dir}loci_fastGear_out \
        -g {input.loci} \
        -o {fastGear_dir}plot_pangenome/pan_fastgear_plot \
        -s True -f pdf -p {input.tree} -z False -y 100 -x 100
        """


 





        



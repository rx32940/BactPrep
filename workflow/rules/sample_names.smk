rule get_sample_names:
    input:
        get_sample_dir
    output:
        out_dir + "all_samples_id.txt"
    script:
        "../scripts/get_sample_names.py"
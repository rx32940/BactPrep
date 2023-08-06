rule get_sample_names:
    input:
        get_sample_dir
    output:
        os.path.join(out_dir , "all_samples_id.txt")
    script:
        "{WORKFLOW}scripts/get_sample_names.py"
rule collect_metadata:
    input:
        sample_table = "config/sample_table.csv",
        wolbachia = "../../resources/wolbachia.xlsx"
    output:
        "results/meta/metadata.csv"
    script:
        "../scripts/collect_metadata.R"

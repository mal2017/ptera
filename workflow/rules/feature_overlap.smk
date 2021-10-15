rule extract_polymorphic_insertions_0:
    input:
        "resources/Tidal_Fly_v1Archive20150930.zip"
    output:
        dgn = temp("results/overlaps/polymorphic_insertions/DGN_flies.zip"),
        labstrain = temp("results/overlaps/polymorphic_insertions/LabStrain_flies.zip"),
        dgrp = temp("results/overlaps/polymorphic_insertions/DGRP_flies.zip"),
        lines = temp("results/overlaps/polymorphic_insertions/CellLines.zip"),
        pools = temp("results/overlaps/polymorphic_insertions/Pool_Flies.zip"),
    shell:
        """
        unzip {input} -d results/overlaps/polymorphic_insertions/
        """

TIDAL_STANDARD_OUTPUTS = ["_Depletion_Annotated.bed","_Depletion_Annotated_TEonly.bed",
"_Depletion_Annotated_TEonly.txt","_Depletion_Annotated_TEonly.txt","_fixed_bin.txt",
"_Inserts_Annotated.bed","_Inserts_Annotated.txt","_map_insertion_depletion.txt",
"_ReadDepletion.txt","_ReadInsertion.txt","_summary.txt","_TE_Indel_genome_plot.pdf"]


checkpoint extract_polymorphic_insertions:
    input:
        "results/overlaps/polymorphic_insertions/{tidal_group}.zip",
    output:
        dir = directory("results/overlaps/polymorphic_insertions/dl/{tidal_group}/")
    shell:
        """
        unzip {input} -d {output}
        """

rule get_reference_copies:
    input:
        "resources/repmasker_dm6_track.txt"
    output:
        "results/overlaps/reference_tes/reference_tes.bed"
    shell:
        """
        cut -f 6,7,8,10,11 {input} | tail -n+2 | awk -F'\t' 'BEGIN {{OFS = FS}} {{t=$4; $4=$5 FS "0" FS t;print}}' | cut -f 1,2,3,4,5,6 > {output}
        """

rule get_per_strain_insertions:
    input:
        ref_copies = rules.bedops_parse_repeatmasker.output,
        poly_dir = rules.extract_polymorphic_insertions.output
    params:
        poly_ins_tsv = rules.extract_polymorphic_insertions.output.dir + "/{tidal_strain}_result/{tidal_strain}_Inserts_Annotated.txt",
        poly_deps_tsv = rules.extract_polymorphic_insertions.output.dir + "/{tidal_strain}_result/{tidal_strain}_Depletion_Annotated.txt",
    output:
        tsv="results/overlaps/polymorphic_insertions/{tidal_group}/{tidal_strain}.tsv"
    script:
        '../scripts/per_strain_ins.R'

def aggregate_known_insertions_within_group(wildcards):
    checkpoint_output = checkpoints.extract_polymorphic_insertions.get(**wildcards).output.dir
    print(checkpoint_output)
    wc_path = os.path.join(checkpoint_output, "{d}_result/")
    print(wc_path)
    x = glob_wildcards(wc_path)
    print(x)
    return expand("results/overlaps/polymorphic_insertions/{tidal_group}/{i}.tsv",
           tidal_group=wildcards.tidal_group,
           i=x.d)

localrules: collect_all_tidal_insertions, collect_known_insertions_within_group

rule collect_known_insertions_within_group:
    input:
        aggregate_known_insertions_within_group
    output:
        tsv="results/overlaps/collected_insertions/{tidal_group}.tsv"
    resources:
        time=60,
        mem=24000,
        cpus=1
    shell:
        "xsv cat rows -d '\t' {input} | tr ',' '\t' > {output}"

rule collect_all_tidal_insertions:
    input:
        expand("results/overlaps/collected_insertions/{g}.tsv",g=["DGN_flies","DGRP_flies","CellLines","Pool_Flies","LabStrain_flies"])
    output:
        "results/overlaps/collected_insertions/insertions_by_strain.tsv.gz"
    resources:
        time=60,
        mem=24000,
        cpus=1
    shell:
        "xsv cat rows -d '\t' {input} | tr ',' '\t' | gzip -c > {output}"

rule feature_overlaps:
    input:
        genes = config.get("TRANSCRIPTOME_GTF"),
        ins = rules.collect_all_tidal_insertions.output
    output:
        tsv="results/overlaps/overlaps.tsv.gz"
    script:
        "../scripts/per_strain_te_gene_overlaps.R"

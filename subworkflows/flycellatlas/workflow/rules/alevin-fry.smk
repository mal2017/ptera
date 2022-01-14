rule update_gtf_for_alevinfry:
    """
    Two important things:
    must have exons for TEs.
    Must have pruning mode set in seqlevels() (see script).
    """
    input:
        gtf = refs_wf("results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf"),
    output:
        gtf = temp("results/tmp_af_txome.gtf")
    script:
        "../scripts/update_gtf_for_alevinfry.R"


rule make_splici:
    """
    https://github.com/COMBINE-lab/alevin-fry
    https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/
    """
    input:
        #genome = refs_wf("results/repeatmasker/genome.fasta.masked"),
        genome = "../../" + config.get("COMBINED_TE_GENOME_FA"),
        gtf = rules.update_gtf_for_alevinfry.output.gtf
    output:
        fasta = "results/transcriptome_splici_fl/transcriptome_splici_fl86.fa",
        t2g_3c = "results/transcriptome_splici_fl/transcriptome_splici_fl86_t2g_3col.tsv",
        t2g = "results/transcriptome_splici_fl/transcriptome_splici_fl86_t2g.tsv"
    resources:
        time=10,
        mem=4000,
        cpus=1
    params:
        read_length = 91,
        flank_trim_length = 5,
        odir="results/transcriptome_splici_fl/",
    script:
        "../scripts/make_splici.R"

rule make_splici_idx:
    input:
        fasta = rules.make_splici.output.fasta
    output:
        idx = directory("results/transcriptome_splici_idx")
    threads:
        24
    resources:
        time=480,
        mem=128000,
        cpus=24
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    shell:
        """
        salmon index -t {input.fasta} \
            -i {output.idx} -p {threads}
        """

rule dl_array_express:
    params:
        r1 = lambda wc: pep.get_sample(wc.sample).fq1_uri[int(wc.ss)],
        r2 = lambda wc: pep.get_sample(wc.sample).fq2_uri[int(wc.ss)],
    output:
        r1 = temp("results/fastqs/{sample}_{ss}_r1.fastq.gz"),
        r2 = temp("results/fastqs/{sample}_{ss}_r2.fastq.gz")
    shell:
        """
        wget -O {output.r1} {params.r1} &&
        wget -O {output.r2} {params.r2}
        """

rule alevin_map:
    """
    https://alevin-fry.readthedocs.io/en/latest/getting_started.html#running-the-alevin-fry-pipeline
    library type choice: https://github.com/COMBINE-lab/salmon/discussions/674
    """
    input:
        idx = rules.make_splici_idx.output.idx,
        r1 = lambda wc: expand("results/fastqs/{{sample}}_{ss}_r1.fastq.gz",ss = range(len(pep.get_sample(wc.sample).subsample_name))),
        r2 = lambda wc: expand("results/fastqs/{{sample}}_{ss}_r2.fastq.gz",ss = range(len(pep.get_sample(wc.sample).subsample_name))),
        tg = rules.make_splici.output.t2g
    output:
        dir = directory("results/alevin-fry/alevin/{sample}")
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    threads:
        24
    resources:
        time=480,
        mem=128000,
        cpus=24
    priority:
        2
    shell:
        """
        salmon alevin -i {input.idx} -p {threads} -l IU --chromiumV3 --sketch \
            -1 {input.r1} \
            -2 {input.r2} \
            --tgMap {input.tg} \
            -o {output.dir}
        """

rule alevin_fry_permit:
    input:
        map = rules.alevin_map.output.dir
    output:
        dir = directory("results/alevin-fry/quant/{sample}")
    threads:
        1
    singularity:
        "docker://quay.io/biocontainers/alevin-fry:0.4.3--h7d875b9_0"
    resources:
        time=20,
        mem=4000,
        cpus=1
    priority:
        2
    shell:
        """
        alevin-fry generate-permit-list -d fw -k -i {input.map} -o {output.dir}
        """

rule alevin_fry_collate:
    input:
        pml = rules.alevin_fry_permit.output.dir,
        map = rules.alevin_map.output.dir
    output:
        touch("results/alevin-fry/quant/{sample}.collate.done")
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/alevin-fry:0.4.3--h7d875b9_0"
    resources:
        time=20,
        mem=4000,
        cpus=12
    priority:
        2
    shell:
        """
        alevin-fry collate -t {threads} -i {input.pml} -r {input.map}
        """

rule alevin_fry_quant:
    input:
        tg3c = rules.make_splici.output.t2g_3c,
        quant = rules.alevin_fry_permit.output.dir,
        collation = rules.alevin_fry_collate.output,
    output:
        dir = directory("results/alevin-fry/quant_res/{sample}")
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/alevin-fry:0.4.3--h7d875b9_0"
    resources:
        time=20,
        mem=4000,
        cpus=12
    priority:
        2
    shell:
        """
        alevin-fry quant -t {threads} \
            -i {input.quant} \
            -o {output.dir} \
            --tg-map {input.tg3c} \
            --resolution cr-like \
            --use-mtx
        """

rule make_r_objs:
    input:
        frydir = rules.alevin_fry_quant.output.dir
    output:
        sce = "results/alevin-fry/sce_objs/{sample}.usa.sce.rds",
        seur = "results/alevin-fry/seur_objs/{sample}.usa.seur.rds",
    script:
        "../scripts/import_alevin_fry.R"

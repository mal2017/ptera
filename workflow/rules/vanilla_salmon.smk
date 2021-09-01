localrules: salmon_decoy_get_ids

rule salmon_decoy_mask_exons:
    """
    We generated a partial decoy-aware transcriptome reference for salmon as described
    in the [salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).
    Briefly, we began by extracting exon genomic coordinates and masking them with bedtools
    v2.30.0 maskfasta.
    """
    input:
        gtf = config.get("TRANSCRIPTOME_GTF"),
        genome = config.get("GENOME_FASTA")
    output:
        exons = "results/quantification/vanilla_salmon_tes_transcripts/decoy/exons.bed",
        genome = "results/quantification/vanilla_salmon_tes_transcripts/decoy/genome.fasta",
        masked = "results/quantification/vanilla_salmon_tes_transcripts/decoy/txome.masked.genome.fasta"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    resources:
        time=10,
        mem=5000,
        cpus=1
    shell:
        """
        awk -v OFS='\\t' '{{if ($3=="exon") {{print $1,$4,$5}}}}' {input.gtf} > {output.exons} &&
        gunzip -c {input.genome} > {output.genome} &&
        bedtools maskfasta -fi {output.genome} -bed {output.exons} -fo {output.masked}
        """

rule salmon_decoy_mashmap_align_txome:
    """
    Next, we applied mashmap v2.0 with options "--pi 80 -s 500" to identify non-genic regions with similarity
    to the transcriptome.
    """
    input:
        masked = rules.salmon_decoy_mask_exons.output.masked,
        txpfile = config.get("FULL_TRANSCRIPT_FASTA")
    output:
        mm = "results/quantification/vanilla_salmon_tes_transcripts/decoy/mashmap.out",
        txlike = "results/quantification/vanilla_salmon_tes_transcripts/decoy/genome_found.sorted.bed",
    threads:
        8
    singularity:
        "docker://quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
    log:
        "results/logs/mashmap/transcripts_and_consensus_tes.txt"
    resources:
        time=20,
        mem=15000,
        cpus=8
    shell:
        """
        mashmap -r {input.masked} -q {input.txpfile} -t {threads} --pi 80 -s 500 -o {output.mm} 2> {log} &&
        awk -v OFS='\\t' '{{print $6,$8,$9}}' {output.mm} | sort -k1,1 -k2,2n - > {output.txlike}
        """

rule salmon_decoy_finalize:
    """
    We extracted the sequence of the transcriptome-similar regions with bedtools getfasta.
    Finally, we concatenated the custom TE and transcriptome reference sequences with the decoy
    sequences to generate the decoy-aware transcriptome reference which we provided to salmon index.
    """
    input:
        txlike = rules.salmon_decoy_mashmap_align_txome.output.txlike,
        masked = rules.salmon_decoy_mask_exons.output.masked,
        txpfile = rules.make_transcripts_and_consensus_tes_fasta.output.fa,
    output:
        txpfile = "results/quantification/vanilla_salmon_tes_transcripts/decoy/txpfile.fasta",
        txlike_merged = "results/quantification/vanilla_salmon_tes_transcripts/decoy/genome_found_merged.bed",
        txlike_fa = "results/quantification/vanilla_salmon_tes_transcripts/decoy/genome_found.fasta",
        decoy = "results/quantification/vanilla_salmon_tes_transcripts/decoy/decoy.fasta",
        gentrome = "results/quantification/vanilla_salmon_tes_transcripts/decoy/gentrome.fasta",
        masked_fai = "results/quantification/vanilla_salmon_tes_transcripts/decoy/txome.masked.genome.fasta.fai"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    resources:
        time=10,
        mem=5000,
        cpus=1
    shell:
        """
        gunzip -c {input.txpfile} > {output.txpfile} &&
        bedtools merge -i {input.txlike} > {output.txlike_merged} &&
        bedtools getfasta -fi {input.masked} -bed {output.txlike_merged} -fo {output.txlike_fa} &&
        awk '{{a=$0; getline;split(a, b, ":");  r[b[1]] = r[b[1]]""$0}} END {{ for (k in r) {{ print k"\\n"r[k] }} }}' {output.txlike_fa} > {output.decoy} &&
        cat {output.txpfile} {output.decoy} > {output.gentrome}
        """

rule salmon_decoy_get_ids:
    input:
        decoy = rules.salmon_decoy_finalize.output.decoy
    output:
        ids = "results/quantification/vanilla_salmon_tes_transcripts/decoy/decoy.txt"
    shell:
        """
        grep ">" {input.decoy} | awk '{{print substr($1,2); }}' > {output.ids}
        """

rule salmon_index_transcripts_and_consensus_tes:
    """
    We indexed the custom transcriptome with salmon v1.5.2 index using "-k <insert k here>."
    """
    input:
        fa = rules.salmon_decoy_finalize.output.gentrome if config.get("SALMON_VANILLA_USE_DECOYS") else rules.make_transcripts_and_consensus_tes_fasta.output.fa,
        decoyids = rules.salmon_decoy_get_ids.output.ids if config.get("SALMON_VANILLA_USE_DECOYS") else rules.make_transcripts_and_consensus_tes_fasta.output.dummy_decoy,
    output:
        directory("results/quantification/vanilla_salmon_tes_transcripts/index/")
    params:
        k = config.get("SALMON_K_PARAM"),
        dec = "--decoys results/quantification/vanilla_salmon_tes_transcripts/decoy/decoy.txt" if config.get("SALMON_VANILLA_USE_DECOYS") else ""
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    threads:
        8
    log:
        "results/logs/salmon_index/transcripts_and_consensus_tes.txt"
    resources:
        time=60,
        mem=20000,
        cpus=8
    shell:
        """
        salmon index -t {input.fa} \
            --index {output} \
            -k {params.k} \
            -p {threads} \
            {params.dec} 2> {log}
        """

rule make_salmon_te_aux_target_file:
    """
    TEs were provided as auxiliary targets to salmon to avoid applying bias correction
    to TE expression estimates.
    """
    input:
        tes = config.get("CONSENSUS_TE_FASTA"),
    output:
        "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.aux.txt"
    shell:
        """
        zcat {input.tes} | grep ">" | tr -d ">" > {output}
        """

rule salmon_quant_se_vanilla:
    """
    We ran salmon v1.5.2 quant directly on trimmed fastqs with options
    "<replace with options actually used>".
    """
    input:
        fq1 = rules.fastp_trim_se.output.r1,
        idx = rules.salmon_index_transcripts_and_consensus_tes.output,
        aux = rules.make_salmon_te_aux_target_file.output,
    output:
        sf = "results/quantification/vanilla_salmon_tes_transcripts/quant/{sample}/quant.sf",
        #sam = temp("results/quantification/vanilla_salmon_tes_transcripts/{sample}/alignments/alignments.sam")
    resources:
        time=60,
        mem=20000,
        cpus=8
    params:
        libtype = config.get("SALMON_VANILLA_LIBTYPE"),
        bootstraps = "--numBootstraps {n}".format(n=config.get("SALMON_VANILLA_BOOTSTRAPS")) if config.get("SALMON_VANILLA_BOOTSTRAPS") else "",
        seqbias = "--seqBias" if config.get("SALMON_VANILLA_SEQBIAS") else "",
        gcbias = "--gcBias" if config.get("SALMON_VANILLA_GCBIAS") else "",
        posbias = "--posBias" if config.get("SALMON_VANILLA_POSBIAS") else "",
        auxtargetfile = "--auxTargetFile {x}".format(x="results/references_and_indices/transcripts_and_consensus_tes/transcripts_and_consensus_tes.aux.txt") if config.get("SALMON_VANILLA_USE_AUXTARGETS") else "",
        softclip = "--softclip" if config.get("SALMON_VANILLA_SOFTCLIP") else "",
    threads:
        8
    log:
        "results/logs/salmon_quant_se_vanilla/{sample}.txt"
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    shell:
        """
        salmon quant --index {input.idx} \
            --libType {params.libtype} \
            -r {input.fq1} \
            {params.bootstraps} \
            {params.seqbias} \
            {params.gcbias} \
            {params.posbias} \
            {params.auxtargetfile} \
            {params.softclip} \
            --writeUnmappedNames \
            -d \
            -p {threads} \
            -o $(dirname {output.sf})/ 2> {log}
        """

rule vanilla_salmon_terminus_group:
    """
    To collapse transcripts to quantifiable clusters for which expression could be accurately estimated,
    we applied terminus v0.1.0 group with default parameters.
    """
    input:
        rules.salmon_quant_se_vanilla.output.sf
    output:
        directory("results/quantification/vanilla_salmon_tes_transcripts/terminus/{sample}/"),
    singularity:
        "docker://quay.io/biocontainers/terminus:0.1.0--hd24f7c9_2"
    resources:
        time=10,
        mem=5000,
        cpus=1
    shell:
        "terminus group -d $(dirname {input}) -o $(dirname {output})"

rule vanilla_salmon_terminus_collapse:
    """
    We used terminus collapse to recalculate expression estimates for all samples.
    """
    input:
        salm = expand("results/quantification/vanilla_salmon_tes_transcripts/quant/{s}/quant.sf", s=SAMPLES),
        term = expand("results/quantification/vanilla_salmon_tes_transcripts/terminus/{s}/", s=SAMPLES)
    output:
        touch("results/quantification/vanilla_salmon_tes_transcripts/terminus.done")
    params:
        od = "results/quantification/vanilla_salmon_tes_transcripts/terminus",
        salm = expand("results/quantification/vanilla_salmon_tes_transcripts/quant/{s}", s=SAMPLES),
    singularity:
        "docker://quay.io/biocontainers/terminus:0.1.0--hd24f7c9_2"
    resources:
        time=10,
        mem=5000,
        cpus=1
    shell:
        """
        terminus collapse -d {params.salm} -o {params.od}
        """

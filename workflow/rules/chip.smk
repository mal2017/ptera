rule macs_chipseq_narrow:
    """
    Silently fails by design. When no peaks are found due to shifting model failure,
    we don't want the whole pipeline to die.
    """
    input:
        ip = rules.dna_cram_to_bam_tmp.output.bam,
        wce = lambda wc: expand("results/mapping/dna_bwa_mem2/{s}.markdup.bam",s=get_background_sample(wc.sample))
    output:
        multiext("results/peaks/chip_macs/narrow/{sample}","_model.r","_peaks.xls","_peaks.narrowPeak","_summits.bed")
    params:
        is_pe = lambda wc: "BAMPE" if is_paired_end(wc.sample) else "AUTO", # TODO - make this pull the info from sample table
        broad = "--broad" if False else "",
        gsize = "dm",
        outdir = "results/peaks/chip_macs/narrow/"
    singularity:
        "docker://quay.io/biocontainers/macs2:2.2.7.1--py37h73a75cf_3"
    priority: 3
    #script:
    #    "../scripts/macs3_chipseq.R"
    shell:
        """
        macs2 callpeak -t {input.ip} -c {input.wce} -f {params.is_pe} \
            -n {wildcards.sample} \
            -g {params.gsize} \
            --outdir {params.outdir} {params.broad} || true
        touch {output}
        """

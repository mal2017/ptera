localrules: get_resources_in_mount

rule get_resources_in_mount:
    """
    Singularity has trouble with files upstream of its mount directory.
    """
    input:
        CONSENSUS_TE_FASTA = "../../" + config.get("CONSENSUS_TE_FASTA"),
        GENOME_FASTA = "../../" + config.get("GENOME_FASTA"),
        FULL_TRANSCRIPT_FASTA = "../../" + config.get("FULL_TRANSCRIPT_FASTA"),
        TRANSCRIPTOME_GTF = "../../" + config.get("TRANSCRIPTOME_GTF"),
        COMBINED_TE_GENOME_FA = "../../" + config.get("COMBINED_TE_GENOME_FA"),
        COMBINED_TE_FEATURE_GTF = "../../" + config.get("COMBINED_TE_FEATURE_GTF")
    output:
        CONSENSUS_TE_FASTA = "results/" + config.get("CONSENSUS_TE_FASTA"),
        GENOME_FASTA = "results/" + config.get("GENOME_FASTA"),
        FULL_TRANSCRIPT_FASTA = "results/" + config.get("FULL_TRANSCRIPT_FASTA"),
        TRANSCRIPTOME_GTF= "results/" + config.get("TRANSCRIPTOME_GTF"),
        COMBINED_TE_GENOME_FA = "results/" + config.get("COMBINED_TE_GENOME_FA"),
        COMBINED_TE_FEATURE_GTF = "results/" + config.get("COMBINED_TE_FEATURE_GTF")
    shell:
        """
        cp {input.CONSENSUS_TE_FASTA} {output.CONSENSUS_TE_FASTA}
        cp {input.GENOME_FASTA} {output.GENOME_FASTA}
        cp {input.FULL_TRANSCRIPT_FASTA} {output.FULL_TRANSCRIPT_FASTA}
        cp {input.TRANSCRIPTOME_GTF} {output.TRANSCRIPTOME_GTF}
        cp {input.COMBINED_TE_GENOME_FA} {output.COMBINED_TE_GENOME_FA}
        cp {input.COMBINED_TE_FEATURE_GTF} {output.COMBINED_TE_FEATURE_GTF}
        """

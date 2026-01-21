import yaml
import os

configfile: "config.yaml"

wildcard_constraints:
    sample="[^/.]+",
    label=r"[0-9]+(?:\.[0-9]+)?"

def get_account_for_jobs(wildcards):
    """Simple alternation based on haplotype"""
    if int(wildcards.hap) == 1:
        return "publicgrp"
    else:
        return "publicgrp"

def seeddotfrac_from_label(label: str, seed: int = 100) -> str:
    frac = float(label) / 100.0
    frac_digits = f"{frac:.10f}".split(".")[1].rstrip("0")
    if not frac_digits:
        frac_digits = "0"
    return f"{seed}.{frac_digits}"

# Build FRAC_LABELS from config
if config.get("dilution", {}).get("enabled", False):
    FRAC_LABELS = [str(x) for x in config["dilution"]["percentages"]]
else:
    FRAC_LABELS = ["100"]

# Build sample list from config
SAMPLES = []
for sample_name, sample_data in config["samples"].items():
    SAMPLES.append(dict(
        sample=sample_name,
        hifi_bam=sample_data["hifi_bam"],
        cifi_bam=sample_data["cifi_bam"],
        enzyme=sample_data["enzyme"]
    ))

def samples_list():
    return [s["sample"] for s in SAMPLES]

def get_sample_data(sample_name):
    return next(s for s in SAMPLES if s["sample"] == sample_name)

def get_enzyme(wildcards):
    return get_sample_data(wildcards.sample)["enzyme"]

# Script paths (relative to workflow)
CIFI2PE_SCRIPT = "scripts/cifi2pe_full_length.py"
CALN50_JS = "scripts/calN50.js"

# External tool paths from config
SINGULARITY_CACHE = config["tools"]["singularity_cache"]

# ---------------- Targets ----------------
rule all:
    input:
        expand("cifi/{sample}.{label}.bam", sample=samples_list(), label=FRAC_LABELS),
        expand("stats/{sample}/{label}/summary.tsv", sample=samples_list(), label=FRAC_LABELS),
        expand("stats/{sample}/{label}/yahs_summary.tsv", sample=samples_list(), label=FRAC_LABELS),
        expand("qc_porec/{sample}/{label}/hap{hap}/bams/{sample}.{label}.hap{hap}.cs.bam",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),

# ---------------- Core steps ----------------

rule hifi_fasta:
    """HiFi BAM -> FASTA (once per sample)"""
    priority: 100
    input:
        bam=lambda w: next(s for s in SAMPLES if s["sample"] == w.sample)["hifi_bam"]
    output:
        fa="hifi/{sample}.hifi.fa"
    threads: 8
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "samtools fasta -@ {threads} {input.bam} > {output.fa}"


rule downsample_cifi_bam:
    """Create label-specific Ci-Fi BAM for pore-c (100% = symlink)."""
    priority: 350
    input:
        bam=lambda w: next(s for s in SAMPLES if s["sample"] == w.sample)["cifi_bam"]
    output:
        bam="cifi/{sample}.{label}.bam",
        bai="cifi/{sample}.{label}.bam.bai"
    params:
        sarg=lambda w: seeddotfrac_from_label(w.label, 100)
    threads: 4
    resources:
        mem_mb=16000, runtime=2 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        r'''
        set -euo pipefail
        if [ "{wildcards.label}" = "100" ]; then
            ln -sf {input.bam} {output.bam}
        else
            samtools view -@ {threads} -b -s {params.sarg} -o {output.bam} {input.bam}
        fi
        samtools index -@ {threads} {output.bam}
        '''


rule cifi_fastq_from_downsampled_bam:
    """Extract FASTQ from downsampled CiFi BAM"""
    priority: 200
    input:
        bam="cifi/{sample}.{label}.bam"
    output:
        fq="cifi/{sample}.{label}.fastq"
    threads: 4
    resources:
        mem_mb=4*1024, runtime=4 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "samtools collate -O -u -@ {threads} {input.bam} | "
        "samtools fastq   -@ {threads} - > {output.fq}"


rule cifi2pe_split:
    """CiFi single FASTQ -> HiC-like PE (R1/R2) with restriction enzyme"""
    priority: 400
    input:  "cifi/{sample}.{label}.fastq"
    output:
        r1="cifi2pe/{sample}.{label}_HiC_R1.fastq",
        r2="cifi2pe/{sample}.{label}_HiC_R2.fastq"
    params:
        out="cifi2pe/{sample}.{label}",
        cutter=get_enzyme
    threads: 1
    resources:
        mem_mb=16000, runtime=12 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "python {CIFI2PE_SCRIPT} --out {params.out} {input} {params.cutter}"

rule hifiasm_dual_scaf:
    """Assemble with hifiasm --dual-scaf (produces hap1/2 ctg GFAs)"""
    priority: 500
    input:
        r1="cifi2pe/{sample}.{label}_HiC_R1.fastq",
        r2="cifi2pe/{sample}.{label}_HiC_R2.fastq",
        hifi="hifi/{sample}.hifi.fa"
    output:
        hap1_gfa="asm/{sample}/{label}/{sample}.{label}.asm.hic.hap1.p_ctg.gfa",
        hap2_gfa="asm/{sample}/{label}/{sample}.{label}.asm.hic.hap2.p_ctg.gfa"
    params:
        pref="asm/{sample}/{label}/{sample}.{label}.asm"
    threads: 32
    resources:
        mem_mb=300*1024, runtime=24 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "hifiasm --dual-scaf -t {threads} -o {params.pref} "
        "--h1 {input.r1} --h2 {input.r2} {input.hifi}"

# ---------------- QC: GFA -> FASTA -> calN50 -> summary ----------------

rule gfa2fa:
    """Convert GFA contigs to FASTA (hap1 or hap2)"""
    priority: 600
    input:
        gfa="asm/{sample}/{label}/{sample}.{label}.asm.hic.hap{hap}.p_ctg.gfa"
    output:
        fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    threads: 4
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "gfatools gfa2fa {input.gfa} > {output.fa}"

rule caln50:
    """Run calN50.js (k8) on each hap FASTA"""
    priority: 700
    input:
        fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        n50="stats/{sample}/{label}/hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"

rule summarize_assembly:
    """
    Parse k8 calN50.js output from hap1/2 and write a tidy TSV with:
    sample, fraction, hap, GS, SZ, NN, N50, L50, AU
    """
    priority: 800
    input:
        hap1="stats/{sample}/{label}/hap1.n50.txt",
        hap2="stats/{sample}/{label}/hap2.n50.txt"
    output:
        tsv="stats/{sample}/{label}/summary.tsv"
    threads: 1
    resources:
        mem_mb=2000, runtime=10, slurm_partition="low", slurm_account="publicgrp"
    shell:
        r'''
        mkdir -p $(dirname {output.tsv})
        echo -e "sample\tfraction\thap\tGS\tSZ\tNN\tN50\tL50\tAU" > {output.tsv}

        parse_one () {{
          IN="$1"; HAP="$2"; SAMP="{wildcards.sample}"; FRAC="{wildcards.label}";
          GS=$(awk -F'\t' '$1=="GS"{{print $2}}' "$IN")
          SZ=$(awk -F'\t' '$1=="SZ"{{print $2}}' "$IN")
          NN=$(awk -F'\t' '$1=="NN"{{print $2}}' "$IN")
          N50=$(awk -F'\t' '$1=="NL" && $2==50{{print $3}}' "$IN")
          L50=$(awk -F'\t' '$1=="NL" && $2==50{{print $4}}' "$IN")
          AU=$(awk -F'\t' '$1=="AU"{{print $2}}' "$IN")
          echo -e "${{SAMP}}\t${{FRAC}}\t${{HAP}}\t${{GS}}\t${{SZ}}\t${{NN}}\t${{N50}}\t${{L50}}\t${{AU}}"
        }}

        parse_one {input.hap1} hap1 >> {output.tsv}
        parse_one {input.hap2} hap2 >> {output.tsv}
        '''


rule porec_nextflow:
    """Run epi2me-labs/wf-pore-c nextflow pipeline for Hi-C contact mapping"""
    priority: 650
    input:
        cifi_bam="cifi/{sample}.{label}.bam",
        ref_fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        bed="porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed",
        pairs="porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.pairs.gz",
        mcool="porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.mcool",
        hic="porec/{sample}/{label}/hap{hap}/hi-c/{sample}.{label}.hap{hap}.hic",
        report="porec/{sample}/{label}/hap{hap}/wf-pore-c-report.html"
    params:
        outdir="porec/{sample}/{label}/hap{hap}",
        sample_alias="{sample}.{label}.hap{hap}",
        cutter=get_enzyme,
        nxf_cache=SINGULARITY_CACHE
    log:
        nf="porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.nextflow.log",
        time="porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.time.txt"
    threads: 16
    resources:
        mem_mb= 32*1024,
        runtime=48 * 60, slurm_partition="low", slurm_account=lambda w: "publicgrp" if int(w.hap) == 1 else "publicgrp"
    shell:
        """
        export NXF_SINGULARITY_CACHEDIR="{params.nxf_cache}"
        export NXF_OPTS='-Xms2g -Xmx4g'
        export NXF_OFFLINE='true'
        export TMPDIR="/scratch/tmp/{wildcards.sample}.{wildcards.label}.hap{wildcards.hap}"
        mkdir -p $TMPDIR
        export NXF_TMP="$TMPDIR"
        export NXF_TEMP="$TMPDIR"
        export NXF_EXECUTOR='local'

        # Create and cd to output directory to isolate nextflow run
        mkdir -p {params.outdir}
        cd {params.outdir}
        echo "Running nextflow in $(pwd)"
        # test relative path to input files
        echo "Input BAM: ../../../../{input.cifi_bam}"
        echo "Input REF: ../../../../{input.ref_fa}"

        /usr/bin/time -v \
            nextflow run epi2me-labs/wf-pore-c \
                -r v1.3.0 \
                -profile singularity \
                --bam ../../../../{input.cifi_bam} \
                --ref ../../../../{input.ref_fa} \
                --cutter {params.cutter} \
                --out_dir . \
                --threads {threads} \
                --minimap2_settings '-ax map-hifi' \
                --paired_end \
                --bed --pairs --hi_c --mcool --coverage \
                --sample "{params.sample_alias}" \
                -with-report   pipeline_report.html \
                -with-timeline pipeline_timeline.html \
                -with-trace    pipeline_trace.txt \
                -with-dag      pipeline_dag.svg \
                1> ../../../../{log.nf} 2> ../../../../{log.time}
        """

rule index_fa:
    """Index FASTA files for yahs"""
    priority: 625
    input:
        fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        fai="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa.fai"
    threads: 1
    resources:
        mem_mb=8 * 1024, runtime=60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "samtools faidx {input.fa}"


rule yahs_scaffold:
    """Scaffold haplotype assemblies with yahs"""
    priority: 750
    input:
        asm_fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa",
        asm_fai="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa.fai",
        porec_bed="porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed"
    output:
        scaffolds="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa",
        agp="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.agp",
        bin="yahs/{sample}/{label}/{sample}.{label}.hap{hap}.bin"
    params:
        prefix="yahs/{sample}/{label}/{sample}.{label}.hap{hap}"
    log:
        "yahs/{sample}/{label}/logs/{sample}.{label}.hap{hap}.yahs.log"
    threads: 32
    resources:
        mem_mb=250*1024, runtime=48 * 60, slurm_partition="low", slurm_account=get_account_for_jobs
    shell:
        """
        yahs \
            -q 0 \
            --no-contig-ec \
            -o {params.prefix} \
            -v 1 \
            {input.asm_fa} \
            {input.porec_bed} \
            2>&1 | tee {log}
        """

rule yahs_caln50:
    """Run calN50.js (k8) on each yahs scaffold FASTA"""
    priority: 850
    input:
        fa="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        n50="stats/{sample}/{label}/yahs_hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition="low", slurm_account="publicgrp"
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"


rule summarize_yahs:
    """
    Parse k8 calN50.js output from yahs scaffolds hap1/2 and write a tidy TSV with:
    sample, fraction, hap, GS, SZ, NN, N50, L50, AU
    """
    priority: 900
    input:
        hap1="stats/{sample}/{label}/yahs_hap1.n50.txt",
        hap2="stats/{sample}/{label}/yahs_hap2.n50.txt"
    output:
        tsv="stats/{sample}/{label}/yahs_summary.tsv"
    threads: 1
    resources:
        mem_mb=2000, runtime=10, slurm_partition="low", slurm_account="publicgrp"
    shell:
        r'''
        mkdir -p $(dirname {output.tsv})
        echo -e "sample\tfraction\thap\tGS\tSZ\tNN\tN50\tL50\tAU" > {output.tsv}

        parse_one () {{
          IN="$1"; HAP="$2"; SAMP="{wildcards.sample}"; FRAC="{wildcards.label}";
          GS=$(awk -F'\t' '$1=="GS"{{print $2}}' "$IN")
          SZ=$(awk -F'\t' '$1=="SZ"{{print $2}}' "$IN")
          NN=$(awk -F'\t' '$1=="NN"{{print $2}}' "$IN")
          N50=$(awk -F'\t' '$1=="NL" && $2==50{{print $3}}' "$IN")
          L50=$(awk -F'\t' '$1=="NL" && $2==50{{print $4}}' "$IN")
          AU=$(awk -F'\t' '$1=="AU"{{print $2}}' "$IN")
          echo -e "${{SAMP}}\t${{FRAC}}\t${{HAP}}\t${{GS}}\t${{SZ}}\t${{NN}}\t${{N50}}\t${{L50}}\t${{AU}}"
        }}

        parse_one {input.hap1} hap1 >> {output.tsv}
        parse_one {input.hap2} hap2 >> {output.tsv}
        '''

rule yahs_index_scaffolds_fa:
    """Index YAHS scaffold FASTA so we can derive chrom.sizes."""
    priority: 740
    input:
        fa="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        fai="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa.fai"
    threads: 1
    resources:
        mem_mb=8*1024, runtime=30, slurm_partition="low", slurm_account=get_account_for_jobs
    shell:
        "samtools faidx {input.fa}"

rule qc_porec_nextflow:
    """Run epi2me-labs/wf-pore-c nextflow pipeline for Hi-C contact mapping on scaffolds"""
    priority: 650
    input:
        cifi_bam="cifi/{sample}.{label}.bam",
        ref_fa="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        bed="qc_porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed",
        pairs="qc_porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.pairs.gz",
        mcool="qc_porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.mcool",
        hic="qc_porec/{sample}/{label}/hap{hap}/hi-c/{sample}.{label}.hap{hap}.hic",
        report="qc_porec/{sample}/{label}/hap{hap}/wf-pore-c-report.html",
        bam="qc_porec/{sample}/{label}/hap{hap}/bams/{sample}.{label}.hap{hap}.cs.bam",
    params:
        outdir="qc_porec/{sample}/{label}/hap{hap}",
        sample_alias="{sample}.{label}.hap{hap}",
        cutter=get_enzyme,
        nxf_cache=SINGULARITY_CACHE
    log:
        nf="qc_porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.nextflow.log",
        time="qc_porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.time.txt"
    threads: 16
    resources:
        mem_mb= 32*1024,
        runtime=48 * 60, slurm_partition="low", slurm_account="publicgrp"
    shell:
        """
        export NXF_SINGULARITY_CACHEDIR="{params.nxf_cache}"
        export NXF_OPTS='-Xms2g -Xmx4g'
        export NXF_OFFLINE='true'
        export TMPDIR="/scratch/tmp/{wildcards.sample}.{wildcards.label}.hap{wildcards.hap}"
        mkdir -p $TMPDIR
        export NXF_TMP="$TMPDIR"
        export NXF_TEMP="$TMPDIR"
        export NXF_EXECUTOR='local'

        # Create and cd to output directory to isolate nextflow run
        mkdir -p {params.outdir}
        cd {params.outdir}
        echo "Running nextflow in $(pwd)"
        # test relative path to input files
        echo "Input BAM: ../../../../{input.cifi_bam}"
        echo "Input REF: ../../../../{input.ref_fa}"

        /usr/bin/time -v \
            nextflow run epi2me-labs/wf-pore-c \
                -r v1.3.0 \
                -profile singularity \
                --bam ../../../../{input.cifi_bam} \
                --ref ../../../../{input.ref_fa} \
                --cutter {params.cutter} \
                --out_dir . \
                --threads {threads} \
                --minimap2_settings '-ax map-hifi' \
                --paired_end \
                --bed --pairs --hi_c --mcool --coverage \
                --sample "{params.sample_alias}" \
                -with-report   pipeline_report.html \
                -with-timeline pipeline_timeline.html \
                -with-trace    pipeline_trace.txt \
                -with-dag      pipeline_dag.svg \
                1> ../../../../{log.nf} 2> ../../../../{log.time}
        """

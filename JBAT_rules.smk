# ============================================================================
# JBAT (Juicebox Assembly Tools) File Preparation Rules
# ============================================================================
# Append these rules to your existing Snakefile
#
# Prerequisites:
#   git clone https://github.com/aidenlab/3d-dna.git
#   wget https://github.com/aidenlab/JuicerTools/releases/download/v3.0.0/juicer_tools.jar
#
# Add to config.yaml:
#   tools:
#     threed_dna: "/path/to/3d-dna"
#     juicer_tools_jar: "/path/to/juicer_tools.jar"
# ============================================================================

# Tool paths from config
THREED_DNA = config["tools"]["threed_dna"]
JUICER_TOOLS_JAR = config["tools"]["juicer_tools_jar"]

# ----------------------------------------------------------------------------
# JBAT file preparation
# ----------------------------------------------------------------------------

rule pairs_to_mnd:
    """
    Convert .pairs.gz from wf-pore-c to merged_nodups.txt format for 3D-DNA.
    
    merged_nodups format (simplified for assembly):
    str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 seq1 mapq2 cigar2 seq2
    """
    input:
        pairs="qc_porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.pairs.gz"
    output:
        mnd="jbat/{sample}/{label}/hap{hap}/merged_nodups.txt"
    threads: 4
    resources:
        mem_mb=32*1024, runtime=4*60, slurm_partition="high", slurm_account="genome-center-grp"
    shell:
        r'''
        mkdir -p $(dirname {output.mnd})
        # Convert 4DN pairs format to merged_nodups.txt format
        # pairs: readID chr1 pos1 chr2 pos2 strand1 strand2
        # mnd:   str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq ...
        zcat {input.pairs} | grep -v "^#" | awk 'BEGIN{{OFS="\t"}} {{
            str1 = ($6 == "+") ? 0 : 16
            str2 = ($7 == "+") ? 0 : 16
            print str1, $2, $3, 0, str2, $4, $5, 1, 60, "-", "-", 60, "-", "-"
        }}' > {output.mnd}
        '''


rule generate_assembly_file:
    """
    Generate .assembly file from yahs scaffold FASTA using 3D-DNA utility.
    """
    input:
        fa="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        assembly="jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.assembly"
    params:
        awk_script=f"{THREED_DNA}/utils/generate-assembly-file-from-fasta.awk"
    threads: 1
    resources:
        mem_mb=8*1024, runtime=30, slurm_partition="high", slurm_account="genome-center-grp"
    shell:
        r'''
        mkdir -p $(dirname {output.assembly})
        awk -f {params.awk_script} {input.fa} > {output.assembly}
        '''


rule jbat_hic:
    """
    Generate JBAT-compatible .hic file using 3D-DNA visualizer.
    This creates the 'assembly' pseudo-chromosome format required by JBAT.
    """
    input:
        assembly="jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.assembly",
        mnd="jbat/{sample}/{label}/hap{hap}/merged_nodups.txt"
    output:
        hic="jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic"
    params:
        threed_dna=THREED_DNA,
        juicer_jar=JUICER_TOOLS_JAR,
        workdir="jbat/{sample}/{label}/hap{hap}",
        prefix="{sample}.{label}.hap{hap}"
    log:
        "jbat/{sample}/{label}/hap{hap}/logs/visualize.log"
    threads: 16
    resources:
        mem_mb=128*1024, runtime=12*60, slurm_partition="high", slurm_account="genome-center-grp"
    shell:
        r'''
        mkdir -p {params.workdir}/logs
        cd {params.workdir}
        
        # Set juicer_tools path for 3D-DNA scripts
        export JUICER_TOOLS_JAR="{params.juicer_jar}"
        
        # Run 3D-DNA visualizer
        bash {params.threed_dna}/visualize/run-assembly-visualizer.sh \
            -p true \
            {params.prefix}.assembly \
            merged_nodups.txt \
            2>&1 | tee logs/visualize.log
        '''


# ----------------------------------------------------------------------------
# Post-curation: Convert reviewed assembly back to FASTA
# ----------------------------------------------------------------------------

rule jbat_post_review:
    """
    After manual curation in JBAT, convert reviewed .assembly to FASTA.
    
    Usage:
      1. Open jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic in Juicebox
      2. Load the .assembly file via Assembly -> Import Map Assembly
      3. Make corrections, then Export Assembly as {sample}.{label}.hap{hap}.review.assembly
      4. Run: snakemake jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.FINAL.fa
    """
    input:
        review="jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.review.assembly",
        orig_fa="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa",
        mnd="jbat/{sample}/{label}/hap{hap}/merged_nodups.txt"
    output:
        final_fa="jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.FINAL.fa"
    params:
        threed_dna=THREED_DNA,
        juicer_jar=JUICER_TOOLS_JAR,
        workdir="jbat/{sample}/{label}/hap{hap}"
    log:
        "jbat/{sample}/{label}/hap{hap}/logs/post_review.log"
    threads: 16
    resources:
        mem_mb=128*1024, runtime=12*60, slurm_partition="high", slurm_account="genome-center-grp"
    shell:
        r'''
        cd {params.workdir}
        export JUICER_TOOLS_JAR="{params.juicer_jar}"
        
        bash {params.threed_dna}/run-asm-pipeline-post-review.sh \
            -r {wildcards.sample}.{wildcards.label}.hap{wildcards.hap}.review.assembly \
            ../../../../{input.orig_fa} \
            merged_nodups.txt \
            2>&1 | tee logs/post_review.log
        
        # 3D-DNA outputs .FINAL.fasta, rename to .FINAL.fa for consistency
        if [ -f *.FINAL.fasta ]; then
            mv *.FINAL.fasta {wildcards.sample}.{wildcards.label}.hap{wildcards.hap}.FINAL.fa
        fi
        '''


# ----------------------------------------------------------------------------
# Target rules
# ----------------------------------------------------------------------------

rule jbat_all:
    """Generate all JBAT files for manual curation"""
    input:
        expand("jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.hic",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),
        expand("jbat/{sample}/{label}/hap{hap}/{sample}.{label}.hap{hap}.assembly",
               sample=samples_list(), label=FRAC_LABELS, hap=[1, 2]),

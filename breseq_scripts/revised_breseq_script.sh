
#!/bin/bash

# === CONFIGURATION === #
RAWREADS="raw_reads"
OUTPUT="processedReads2"
SAMPLES_FILE="samples.txt"
SUBSAMPLED="${OUTPUT}/subsampled"
THREADS=32
MINLEN=100
SUBSAMPLE_SIZE=500000

# === DATABASE PATHS === #
CENTRIFUGE_WOL_DB="/galaxy/biodatabases/web_of_life/centrifuge/WoLr1"
CENTRIFUGE_DEFAULT_DB="/galaxy/biodatabases/centrifuge/09262018/p+h+v"
KRAKEN2_DB="/galaxy/biodatabases/kraken2-default"
BRESEQ_REF="../Ecoli_NIST0056.gbk"

# === SETUP === #
mkdir -p "$SUBSAMPLED"
mkdir -p results/{centrifugeWoL,centrifuge,kraken2,breseq}

# === FUNCTIONS === #


generate_sample_list() {
    echo "Generating sample list from $RAWREADS..."
    cd "$RAWREADS" || exit 1

    # Clear existing samples.txt and regenerate
    > "../$SAMPLES_FILE"
    for f in *_R1.fastq.gz; do
        [[ -f "$f" ]] || continue
        echo "${f%%_R1.fastq.gz}" >> "../$SAMPLES_FILE"
    done

    cd - || exit 1
}

trim_reads() {
    echo "Starting trimming..."
    eval "$(conda shell.bash hook)"
    conda activate breseq0_35_7_env

    cd "$RAWREADS" || exit 1
    echo -e "SampleID\tReadCount" > "../$OUTPUT/processedReads2.txt"

    while read -r sampleID; do
        echo "Trimming $sampleID..."

        in1="${sampleID}_R1.fastq.gz"
        in2="${sampleID}_R2.fastq.gz"

        if [[ ! -f "${in1}" || ! -f "${in2}" ]]; then
            echo "Missing input files for $sampleID — skipping."
            continue
        fi

        bbduk.sh -Xmx84g \
            in1="${in1}" \
            in2="${in2}" \
            out1="../$OUTPUT/${sampleID}_R1.fastq.gz" \
            out2="../$OUTPUT/${sampleID}_R2.fastq.gz" \
            t="$THREADS" \
            qtrim=rl trimq=20 minlen="$MINLEN" maxns=0 rieb=t \
            |& tee -a "../$OUTPUT/${sampleID}_bbduk_rpt.txt"

        readCount=$(zgrep "^@" "../$OUTPUT/${sampleID}_R1.fastq.gz" | wc -l)
        echo -e "${sampleID}\t${readCount}" >> "../$OUTPUT/processedReads2.txt"
    done < "../$SAMPLES_FILE"

    cd - || exit 1
}

subsample_reads() {
    echo "Starting subsampling..."
    cd "$OUTPUT" || exit 1

    while read -r sampleID; do
        echo "Subsampling $sampleID..."
        zcat "${sampleID}_R1.fastq.gz" | seqtk sample -s100 - "$SUBSAMPLE_SIZE" | gzip > "$SUBSAMPLED/${sampleID}_500K_R1.fastq.gz"
        zcat "${sampleID}_R2.fastq.gz" | seqtk sample -s100 - "$SUBSAMPLE_SIZE" | gzip > "$SUBSAMPLED/${sampleID}_500K_R2.fastq.gz"
        zcat "$SUBSAMPLED/${sampleID}_500K_R1.fastq.gz" "$SUBSAMPLED/${sampleID}_500K_R2.fastq.gz" > "$SUBSAMPLED/${sampleID}_combined_1M.fastq"
    done < "../$SAMPLES_FILE"

    cd - || exit 1
}

run_centrifuge_wol() {
    echo "Running Centrifuge (WoL)..."
    while read -r sampleID; do
        centrifuge -x "$CENTRIFUGE_WOL_DB" \
            -1 <(zcat "$OUTPUT/${sampleID}_R1.fastq.gz") \
            -2 <(zcat "$OUTPUT/${sampleID}_R2.fastq.gz") \
            -p "$THREADS" \
            -S /dev/null \
            --report-file "results/centrifugeWoL/${sampleID}_centrifuge.txt"

        (head -n 1 "results/centrifugeWoL/${sampleID}_centrifuge.txt" &&
         tail -n +2 "results/centrifugeWoL/${sampleID}_centrifuge.txt" | sort -nr -k5 -t$'\t') \
         > "results/centrifugeWoL/${sampleID}_centrifuge_sorted.txt"
    done < "$SAMPLES_FILE"
}

run_centrifuge_default() {
    echo "Running Centrifuge (Default DB)..."
    while read -r sampleID; do
        centrifuge -x "$CENTRIFUGE_DEFAULT_DB" \
            -1 <(zcat "$OUTPUT/${sampleID}_R1.fastq.gz") \
            -2 <(zcat "$OUTPUT/${sampleID}_R2.fastq.gz") \
            -p "$THREADS" \
            -S /dev/null \
            --report-file "results/centrifuge/${sampleID}_centrifuge.txt"

        (head -n 1 "results/centrifuge/${sampleID}_centrifuge.txt" &&
         tail -n +2 "results/centrifuge/${sampleID}_centrifuge.txt" | sort -nr -k5 -t$'\t') \
         > "results/centrifuge/${sampleID}_centrifuge_sorted.txt"
    done < "$SAMPLES_FILE"
}

run_kraken2() {
    echo "Running Kraken2..."
    export PATH="/environments/builds/kraken2/scripts:$PATH"

    while read -r sampleID; do
        kraken2 --db "$KRAKEN2_DB" \
            --threads "$THREADS" \
            --use-names \
            --confidence 0.05 \
            --report "results/kraken2/${sampleID}_k2_report.txt" \
            --paired \
            "$OUTPUT/${sampleID}_R1.fastq.gz" "$OUTPUT/${sampleID}_R2.fastq.gz" \
            --output "results/kraken2/${sampleID}.txt"

        cut -f1-3,6-8 "results/kraken2/${sampleID}_k2_report.txt" > "results/kraken2/${sampleID}_stdK2_report.txt"
    done < "$SAMPLES_FILE"
}

run_breseq() {
    echo "Running breseq..."
    eval "$(conda shell.bash hook)"
    conda activate breseq0_35_7_env

    cd "$OUTPUT" || exit 1
    while read -r sampleID; do
        R1="${sampleID}_R1.fastq.gz"
        R2="${sampleID}_R2.fastq.gz"

        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "Missing FASTQ files for $sampleID — skipping."
            continue
        fi

        breseq -r "$BRESEQ_REF" -p -j "$THREADS" "$R1" "$R2" \
            -o "../results/breseq/${sampleID}"
    done < "../$SAMPLES_FILE"
    cd - || exit 1
}


generate_staged_bam() {
    echo "Running staged Bowtie2 alignment for cn.MOPS..."
    eval "$(conda shell.bash hook)"
    conda activate breseq0_35_7_env # Replace with your actual env

    REFERENCE_INDEX="/home/william-ackerman/Desktop/cnMOPS/reference"
    RAW_READS="/home/william-ackerman/Desktop/shortcut-to-nas/Microbiology_Processing/Stabryla_Data/Motility_AMR/raw_reads"
    BAM_OUTPUT="results/bam"
    THREADS=32

    mkdir -p "$BAM_OUTPUT"

    # Build Bowtie2 index if not already built
    #if [ ! -f "${REFERENCE_INDEX}.1.bt2" ]; then
        #bowtie2-build "$BRESEQ_REF" "$REFERENCE_INDEX"
    #fi

    while read -r sampleID; do
        echo "Processing $sampleID..."

        OUTDIR="$BAM_OUTPUT/$sampleID"
        mkdir -p "$OUTDIR"

        for READ in R1 R2; do
            # Stage 1
            bowtie2 -t --no-unal -p "$THREADS" -L 31 --ma 1 --mp 3 --np 0 \
                --rdg 2,3 --rfg 2,3 --ignore-quals --local -i S,1,0.25 \
                --score-min L,1,0.9 -k 2000 --reorder \
                -x "$REFERENCE_INDEX" \
                -U "$RAW_READS/${sampleID}_${READ}.fastq.gz" \
                -S "$OUTDIR/${sampleID}_${READ}.stage1.sam" \
                --un "$OUTDIR/${sampleID}_${READ}.stage1.unmatched.fastq"

            # Stage 2
            bowtie2 -t --no-unal -p "$THREADS" -L 19 --ma 1 --mp 3 --np 0 \
                --rdg 2,3 --rfg 2,3 --ignore-quals --local -i S,1,0.25 \
                --score-min L,6,0.2 -k 2000 --reorder \
                -x "$REFERENCE_INDEX" \
                -U "$OUTDIR/${sampleID}_${READ}.stage1.unmatched.fastq" \
                -S "$OUTDIR/${sampleID}_${READ}.stage2.matched.sam"
        done

        # Convert each SAM to BAM
        samtools view -bS "$OUTDIR/${sampleID}_R1.stage1.sam" -o "$OUTDIR/${sampleID}_R1.stage1.bam"
        samtools view -bS "$OUTDIR/${sampleID}_R1.stage2.matched.sam" -o "$OUTDIR/${sampleID}_R1.stage2.bam"
        samtools view -bS "$OUTDIR/${sampleID}_R2.stage1.sam" -o "$OUTDIR/${sampleID}_R2.stage1.bam"
        samtools view -bS "$OUTDIR/${sampleID}_R2.stage2.matched.sam" -o "$OUTDIR/${sampleID}_R2.stage2.bam"

	# Merge BAMs
	samtools merge "$OUTDIR/${sampleID}.merged.bam" \
    		"$OUTDIR/${sampleID}_R1.stage1.bam" \
    		"$OUTDIR/${sampleID}_R1.stage2.bam" \
    		"$OUTDIR/${sampleID}_R2.stage1.bam" \
   		"$OUTDIR/${sampleID}_R2.stage2.bam"

	# Sort and index
	samtools sort -o "$OUTDIR/${sampleID}.sorted.bam" "$OUTDIR/${sampleID}.merged.bam"
	samtools index "$OUTDIR/${sampleID}.sorted.bam"

        # Optional cleanup
        #rm "$OUTDIR"/*.sam "$OUTDIR"/*.unmatched.fastq
    done < "$SAMPLES_FILE"
}

validate_samples() {
    echo "Validating sample list..."
    while read -r sampleID; do
        if [[ ! -f "$OUTPUT/${sampleID}_R1.fastq.gz" || ! -f "$OUTPUT/${sampleID}_R2.fastq.gz" ]]; then
            echo "Warning: Missing files for $sampleID"
        fi
    done < "$SAMPLES_FILE"
}

# THIS NEEDS WORK


# === MAIN PIPELINE === #
main() {
    generate_sample_list
    # trim_reads
    # subsample_reads
    # run_centrifuge_wol
    # run_centrifuge_default
    # run_kraken2
    # validate_samples
    # run_breseq
    generate_staged_bam
}

main


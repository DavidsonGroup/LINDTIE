/*
Module      : assembly
Description : Assembles transcripts for each case sample and merges with reference transcriptome.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

//------------------------------------------------------------------------------
/*
  Process: align_raw_reads_to_hg38
    - Aligns the raw reads to the reference genome and separates the mapped and unmapped reads.
Container:
      - bioconda::minimap2=2.30
      - bioconda::samtools=1.22
*/
process align_raw_reads_to_hg38 {
    
    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output/01-Assembly/", mode: 'copy', pattern: '*.{log,bam,bam.bai}'

    container 'oras://community.wave.seqera.io/library/minimap2_samtools:e98addfcfd60e8e7'

    input:
      tuple val(sample_id), path(reads)
      path hg38_fasta
      path hg38_splice_junctions

    output:
      tuple val(sample_id), path("confident_mapped.bam"), emit: confident_mapped_bam
      tuple val(sample_id), path("confident_mapped.bam.bai"), emit: confident_mapped_bam_bai
      tuple val(sample_id), path("reads_all_sorted.bam"), emit: all_mapped_bam
      tuple val(sample_id), path("reads_all_sorted.bam.bai"), emit: all_mapped_bam_bai
      tuple val(sample_id), path("rescued.fastq"), emit: rescued_fastq
      tuple val(sample_id), path("read_counts_summary.log"), emit: read_counts_summary, optional: true

    script:
    """
    # 1. Align reads to reference genome
    minimap2 -t ${task.cpus} -ax splice --secondary=no --junc-bed ${hg38_splice_junctions} \\
      ${hg38_fasta} ${reads} | \\
      samtools sort -@ ${task.cpus} -o reads_all_sorted.bam

    samtools index reads_all_sorted.bam

    # 2. Two-Pass Routing
    # Pass 1: collect read names that have any supplementary alignment
    samtools view reads_all_sorted.bam | \
      awk 'int(\$2 / 2048) % 2 { print \$1 }' | sort -u > supp_ids.txt

    # Pass 2: route reads; any read with supplementary alignment goes to rescued
    awk -v mode="${params.assembly_mode}" '
      FNR==NR { has_supp[\$1]=1; next }
      BEGIN {
        n_total=0; n_confident=0; n_rescued=0;
        r_unmapped=0; r_supp=0; r_low_mapq=0; r_indel=0; r_clip=0;
      }
      
      /^@/ { 
        print > "confident.sam"
        print > "rescued.sam"
        next 
      }
      
      {
        n_total++
        total_reads[\$1]=1
        flag = \$2
        mapq = \$5
        cigar = \$6
        
        is_unmapped = int(flag / 4) % 2
        is_supp = int(flag / 2048) % 2
        is_low_mapq = (mapq < 20)
        
        has_indel = 0
        is_clipped = 0
        
        if (!is_unmapped) {
            tmp = cigar
            clipped_bases = 0
            while (match(tmp, /[0-9]+[IDNSM]/)) {
                op = substr(tmp, RLENGTH, 1)
                len = substr(tmp, RSTART, RLENGTH - 1) + 0
                if ((op == "I" || op == "D") && len > 15) has_indel = 1
                if (op == "S") clipped_bases += len
                tmp = substr(tmp, RSTART + RLENGTH)
            }
            if (clipped_bases > 50) is_clipped = 1
        }

        if (has_supp[\$1]) {
            print > "rescued.sam"
            n_rescued++
            r_supp++
            rescued_reads[\$1]=1
            rescue_supp[\$1]=1
        } else if (is_unmapped || is_supp || is_low_mapq || has_indel || is_clipped) {
            print > "rescued.sam"
            n_rescued++
            if (is_unmapped) r_unmapped++
            else if (is_supp) r_supp++
            else if (is_low_mapq) r_low_mapq++
            else if (has_indel) r_indel++
            else if (is_clipped) r_clip++
            rescued_reads[\$1]=1
            if (is_unmapped) rescue_unmapped[\$1]=1
            else if (is_supp) rescue_supp[\$1]=1
            else if (is_low_mapq) rescue_low_mapq[\$1]=1
            else if (has_indel) rescue_indel[\$1]=1
            else if (is_clipped) rescue_clip[\$1]=1
        } else {
            print > "confident.sam"
            n_confident++
            confident_reads[\$1]=1
        }
      }
    
    END {
        # Unique read-name counts
        for (r in total_reads) u_total++
        for (r in confident_reads) u_confident++
        for (r in rescued_reads) u_rescued++
        for (r in rescue_unmapped) u_r_unmapped++
        for (r in rescue_supp) u_r_supp++
        for (r in rescue_low_mapq) u_r_low_mapq++
        for (r in rescue_indel) u_r_indel++
        for (r in rescue_clip) u_r_clip++

        # --- LOGGING SECTION ---
        print "Sample Alignment Summary" > "read_counts_summary.log"
        print "Mode: " mode >> "read_counts_summary.log"
        print "Alignment-level counts (SAM records)" >> "read_counts_summary.log"
        print "Total Alignments Processed: " n_total >> "read_counts_summary.log"
        print "--------------------------------" >> "read_counts_summary.log"

        if (n_total > 0) {
            if (mode == "hybrid") {
                # Hybrid: Distinct split between Confident (StringTie) and Rescued (RNABloom)
                print "Sent to StringTie2 (Confident): " n_confident " (" (n_confident/n_total)*100 "%)" >> "read_counts_summary.log"
                print "Sent to RNA-Bloom2 (Rescued):   " n_rescued " (" (n_rescued/n_total)*100 "%)" >> "read_counts_summary.log"
                
                print "--------------------------------" >> "read_counts_summary.log"
                print "Rescue Reasons (Alignment-level):" >> "read_counts_summary.log"
                print "  Unmapped:        " r_unmapped >> "read_counts_summary.log"
                print "  Supplementary:   " r_supp >> "read_counts_summary.log"
                print "  Low MAPQ (<20):  " r_low_mapq >> "read_counts_summary.log"
                print "  Large Indel:     " r_indel >> "read_counts_summary.log"
                print "  Soft Clipped:    " r_clip >> "read_counts_summary.log"
            } 
            else if (mode == "denovo_subset" || mode == "ref_guided") {
                # Denovo Subset / Ref Guided: We use ALL mapped reads (confident + poor quality)
                # We calculate "Mapped" as Total - Unmapped
                n_mapped = n_total - r_unmapped
                print "Sent to StringTie2 (All Mapped): " n_mapped " (" (n_mapped/n_total)*100 "%)" >> "read_counts_summary.log"
                print "Discarded (Unmapped):            " r_unmapped " (" (r_unmapped/n_total)*100 "%)" >> "read_counts_summary.log"
                
                if (mode == "denovo_subset") {
                  print "(Note: These are the 'Remainder' reads not selected for De Novo)" >> "read_counts_summary.log"
                }
            }
            if (mode == "hybrid") {
            print "--------------------------------" >> "read_counts_summary.log"
            print "Read-level counts (unique read names)" >> "read_counts_summary.log"
            print "Unique Read Names (Total): " u_total >> "read_counts_summary.log"
            if (u_total > 0) {
              print "Unique Reads to StringTie2:      " u_confident " (" (u_confident/u_total)*100 "%)" >> "read_counts_summary.log"
              print "Unique Reads to RNA-Bloom2:      " u_rescued " (" (u_rescued/u_total)*100 "%)" >> "read_counts_summary.log"
            }
              print "Unique Rescue Reasons (Read-level):" >> "read_counts_summary.log"
              print "  Unmapped:        " u_r_unmapped >> "read_counts_summary.log"
              print "  Supplementary:   " u_r_supp >> "read_counts_summary.log"
              print "  Low MAPQ (<20):  " u_r_low_mapq >> "read_counts_summary.log"
              print "  Large Indel:     " u_r_indel >> "read_counts_summary.log"
              print "  Soft Clipped:    " u_r_clip >> "read_counts_summary.log"
            }
        } else {
            print "WARNING: No reads found in input." >> "read_counts_summary.log"
        }
      }
    ' supp_ids.txt <(samtools view -h reads_all_sorted.bam)

    # 3. Convert outputs
    samtools view -b -@ ${task.cpus} confident.sam > confident_mapped.bam
    samtools index confident_mapped.bam
    rm confident.sam

    samtools view -b -@ ${task.cpus} rescued.sam > rescued.bam
    samtools fastq -@ ${task.cpus} -n rescued.bam > rescued.fastq
    rm rescued.sam rescued.bam
    """
}

//------------------------------------------------------------------------------
/*
  Process: ref_guided_assembly
    - Use aligned reads (either all or confident) to perform reference-guided assembly using StringTie2.
Container:
        - bioconda::gffread=0.12.7
        - bioconda::stringtie=2.2.3
*/
process ref_guided_assembly {
    
    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output/01-Assembly/", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/gffread_stringtie:12cc4a646c48604f'

    input:
      tuple val(sample_id), path(mapped_bam)
      path tx_annotation
      path hg38_fasta

    output:
      tuple val(sample_id), path("stringtie2_novel_transcripts.gtf"), emit: stringtie2_gtf
      tuple val(sample_id), path("stringtie2_assembly.fa"), emit: stringtie2_assembled_fa

    script:
    """
    # Assemble transcripts
    stringtie ${mapped_bam} -G ${tx_annotation} -p ${task.cpus} -L -o stringtie2.gtf

    # Extract novel transcripts & all exons belonging to those novel transcripts
    awk '\$3=="transcript" && \$0 !~ /reference_id/ {
            print
            tid=""
            for(i=1;i<=NF;i++){
                if(\$i=="transcript_id"){
                    tid=\$(i+1)
                    gsub(/"/,"",tid)
                    gsub(/;/,"",tid)
                    keep[tid]=1
                }
            }
        }
        \$3=="exon" {
            tid=""
            for(i=1;i<=NF;i++){
                if(\$i=="transcript_id"){
                    tid=\$(i+1)
                    gsub(/"/,"",tid)
                    gsub(/;/,"",tid)
                    if(keep[tid]==1) print
                }
            }
        }' stringtie2.gtf > stringtie2_novel_transcripts.gtf

    # Convert GTF to FASTA
    gffread stringtie2_novel_transcripts.gtf -g ${hg38_fasta} -w stringtie2_assembly.fa
    """
}

//------------------------------------------------------------------------------
/*
  Process: subset_reads
    - Splits reads into two sets using a 2-step approach:
      1. Sample the target count to create the subset (reformat.sh).
      2. Filter the original reads against the subset to find the remainder (filterbyname.sh).
  Container:
      - bioconda::bbmap=39.52
*/
process subset_reads {
    tag "${sample_id}"
    label 'process_medium'

    container 'oras://community.wave.seqera.io/library/bbmap:39.52--be327f090dc08662'

    input:
      tuple val(sample_id), path(reads)
      val target_count

    output:
      tuple val(sample_id), path("${sample_id}_subset.fastq"), emit: subset_fastq
      tuple val(sample_id), path("${sample_id}_remainder.fastq"), emit: remainder_fastq

    script:
    """
# -------------------------------------------------------
    # Step 1: Generate the Random Subset
    # -------------------------------------------------------
    # samplereads=${target_count} : exact number of reads to sample
    reformat.sh \\
        in=${reads} \\
        out=${sample_id}_subset.fastq \\
        samplereads=${target_count} \\
        qin=33 \\
        overwrite=true

    # -------------------------------------------------------
    # Step 2: Extract the Remainder
    # -------------------------------------------------------
    # use filterbyname.sh to find reads that are NOT in the subset.
    # names=... : The file containing the names to filter against (the subset)
    # include=f : "False", to get reads that are NOT in the subset file
    
    filterbyname.sh \\
        -Xmx16g \\
        in=${reads} \\
        out=${sample_id}_remainder.fastq \\
        names=${sample_id}_subset.fastq \\
        include=f \\
        qin=33 \\
        overwrite=true
    """
}

//------------------------------------------------------------------------------
/*
  Process: de_novo_assembly
    - Use reads (raw or rescued) to perform de novo assembly using RNA-Bloom2.
Container:
      - bioconda::rnabloom=2.0.1
*/
process de_novo_assembly {

    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/rnabloom:2.0.1--1a308388e7330445'

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("01-Assembly/rnabloom.transcripts.fa"), emit: rnabloom_assembled_fa
      tuple val(sample_id), path("01-Assembly/denovo_read_counts.log"), emit: read_counts_log, optional: true

    script:
    // Automatically removes existing dash (if any) and adds a fresh one
    // This makes inputs 'lrpb' and '-lrpb' both valid
    def bloom_mode = "-" + params.rnabloom2_preset.replaceFirst(/^-/, '')
    """
    ## Create directories
    mkdir -p 01-Assembly

    ## ---------------------------------------------------------
    ## 0. De-duplicate Reads (by read name)
    ## ---------------------------------------------------------
    READS_INPUT="${reads}"
    FIRST_CHAR=\$(head -n 1 "${reads}" | cut -c 1)
    
    if [ "\$FIRST_CHAR" == "@" ]; then
        # FASTQ: remove duplicate read names (keep first)
        awk 'NR%4==1 { h=\$0; sub(/^@/,"",h); dup=seen[h]++ }
             { if (!dup) print }' "${reads}" > dedup_reads.fastq
        READS_INPUT="dedup_reads.fastq"
    elif [ "\$FIRST_CHAR" == ">" ]; then
        # FASTA: remove duplicate read names (keep first)
        awk '/^>/ { h=\$0; sub(/^>/,"",h); dup=seen[h]++ }
             { if (!dup) print }' "${reads}" > dedup_reads.fasta
        READS_INPUT="dedup_reads.fasta"
    fi

    ## ---------------------------------------------------------
    ## 1. Count Reads for Logging
    ## ---------------------------------------------------------
    # Count raw input reads (before dedup)
    RAW_FIRST_CHAR=\$(head -n 1 "${reads}" | cut -c 1)

    if [ "\$RAW_FIRST_CHAR" == "@" ]; then
        RAW_LINE_COUNT=\$(wc -l < "${reads}")
        RAW_READ_COUNT=\$((RAW_LINE_COUNT / 4))
    elif [ "\$RAW_FIRST_CHAR" == ">" ]; then
        RAW_READ_COUNT=\$(grep -c "^>" "${reads}")
    else
        RAW_READ_COUNT="Unknown"
    fi

    # Detect format by checking the first character of the file
    # @ = FASTQ, > = FASTA
    FIRST_CHAR=\$(head -n 1 "\${READS_INPUT}" | cut -c 1)
    
    if [ "\$FIRST_CHAR" == "@" ]; then
        # FASTQ: Count total lines and divide by 4
        LINE_COUNT=\$(wc -l < "\${READS_INPUT}")
        READ_COUNT=\$((LINE_COUNT / 4))
    elif [ "\$FIRST_CHAR" == ">" ]; then
        # FASTA: Count lines starting with '>'
        READ_COUNT=\$(grep -c "^>" "\${READS_INPUT}")
    else
        READ_COUNT="Unknown"
    fi

    # Write to log file
    echo "Sample ID: ${sample_id}" > 01-Assembly/denovo_read_counts.log
    echo "Assembly Mode: ${params.assembly_mode}" >> 01-Assembly/denovo_read_counts.log
    echo "Subset Count: ${params.subset_count}" >> 01-Assembly/denovo_read_counts.log
    echo "--------------------------------" >> 01-Assembly/denovo_read_counts.log
    echo "Reads Input (Raw): \$RAW_READ_COUNT" >> 01-Assembly/denovo_read_counts.log
    echo "Reads Input (Dedup): \$READ_COUNT" >> 01-Assembly/denovo_read_counts.log

    ## ---------------------------------------------------------
    ## 2. Run RNABloom2
    ## ---------------------------------------------------------
    ## Run RNABloom2 - params.rnabloom2_preset is optional only for '-lrpb' for PacBio data.
    rnabloom ${bloom_mode} -long "\${READS_INPUT}" -t ${task.cpus} -outdir 01-Assembly
    """
}

//------------------------------------------------------------------------------
/*
  Process: merge_refTrans_assembly
    - Concatenates the reference transcriptome with the assembled transcripts.
*/
process merge_refTrans_assembly {
    
    tag "${sample_id}"
    label 'process_short'

    input:
      tuple val(sample_id), path(assemblies)
      path trans_fasta

    output:
      // Ref + Novel
      tuple val(sample_id), path("merged_refTrans_assembly.fa"), emit: merged_ref, optional: true
      // Novel Assemblies Only (Combined)
      tuple val(sample_id), path("novel_assemblies_only.fa"), emit: novel_only, optional: true

    script:
    """
    # 1. Concatenate all novel assemblies (handles Hybrid, Denovo, or Ref-Guided)
    cat ${assemblies} > novel_assemblies_only.fa

    # 2. Create the full merged reference
    cat ${trans_fasta} novel_assemblies_only.fa > merged_refTrans_assembly.fa
    """
}
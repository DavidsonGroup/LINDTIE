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

    publishDir "${sample_id}_output/01-Assembly/", mode: 'copy', pattern: '*.{log,bam}'

    container 'oras://community.wave.seqera.io/library/minimap2_samtools:e98addfcfd60e8e7'

    input:
      tuple val(sample_id), path(reads)
      path hg38_fasta

    output:
      tuple val(sample_id), path("confident_mapped.bam"), emit: confident_mapped_bam
      tuple val(sample_id), path("reads_all_sorted.bam"), emit: all_mapped_bam
      tuple val(sample_id), path("rescued.fastq"), emit: rescued_fastq
      tuple val(sample_id), path("read_counts_summary.log"), emit: read_counts_summary, optional: true

    script:
    """
    # 1. Align reads to reference genome
    minimap2 -t ${task.cpus} -ax splice --secondary=no \\
      ${hg38_fasta} ${reads} | \\
      samtools sort -@ ${task.cpus} -o reads_all_sorted.bam

    samtools index reads_all_sorted.bam

    # 2. Single-Pass Routing
    # Pass 'mode' variable to AWK
    samtools view -h reads_all_sorted.bam |
    awk -v mode="${params.assembly_mode}" '
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

        if (is_unmapped || is_supp || is_low_mapq || has_indel || is_clipped) {
            print > "rescued.sam"
            n_rescued++
            if (is_unmapped) r_unmapped++
            else if (is_supp) r_supp++
            else if (is_low_mapq) r_low_mapq++
            else if (has_indel) r_indel++
            else if (is_clipped) r_clip++
        } else {
            print > "confident.sam"
            n_confident++
        }
      }
    
    END {
        # --- LOGGING SECTION ---
        print "Sample Alignment Summary" > "read_counts_summary.log"
        print "Mode: " mode >> "read_counts_summary.log"
        print "Total Reads Processed: " n_total >> "read_counts_summary.log"
        print "--------------------------------" >> "read_counts_summary.log"

        if (n_total > 0) {
            if (mode == "hybrid") {
                # Hybrid: Distinct split between Confident (StringTie) and Rescued (RNABloom)
                print "Sent to StringTie2 (Confident): " n_confident " (" (n_confident/n_total)*100 "%)" >> "read_counts_summary.log"
                print "Sent to RNA-Bloom2 (Rescued):   " n_rescued " (" (n_rescued/n_total)*100 "%)" >> "read_counts_summary.log"
                
                print "--------------------------------" >> "read_counts_summary.log"
                print "Rescue Reasons:" >> "read_counts_summary.log"
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
        } else {
            print "WARNING: No reads found in input." >> "read_counts_summary.log"
        }
      }
    '

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
      path gencode_annotation
      path hg38_fasta

    output:
      tuple val(sample_id), path("stringtie2_novel_transcripts.gtf"), emit: stringtie2_gtf
      tuple val(sample_id), path("stringtie2_assembly.fa"), emit: stringtie2_assembled_fa

    script:
    """
    # Assemble transcripts
    stringtie ${mapped_bam} -G ${gencode_annotation} -p ${task.cpus} -L -o stringtie2.gtf

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
    ## 1. Count Reads for Logging
    ## ---------------------------------------------------------
    # Detect format by checking the first character of the file
    # @ = FASTQ, > = FASTA
    FIRST_CHAR=\$(head -n 1 ${reads} | cut -c 1)
    
    if [ "\$FIRST_CHAR" == "@" ]; then
        # FASTQ: Count total lines and divide by 4
        LINE_COUNT=\$(wc -l < ${reads})
        READ_COUNT=\$((LINE_COUNT / 4))
    elif [ "\$FIRST_CHAR" == ">" ]; then
        # FASTA: Count lines starting with '>'
        READ_COUNT=\$(grep -c "^>" ${reads})
    else
        READ_COUNT="Unknown"
    fi

    # Write to log file
    echo "Sample ID: ${sample_id}" > 01-Assembly/denovo_read_counts.log
    echo "Assembly Mode: ${params.assembly_mode}" >> 01-Assembly/denovo_read_counts.log
    echo "Subset Count: ${params.subset_count}" >> 01-Assembly/denovo_read_counts.log
    echo "--------------------------------" >> 01-Assembly/denovo_read_counts.log
    echo "Reads Input to RNABloom2: \$READ_COUNT" >> 01-Assembly/denovo_read_counts.log

    ## ---------------------------------------------------------
    ## 2. Run RNABloom2
    ## ---------------------------------------------------------
    ## Run RNABloom2 - params.rnabloom2_preset is optional only for '-lrpb' for PacBio data.
    rnabloom ${bloom_mode} -long ${reads} -t ${task.cpus} -outdir 01-Assembly
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
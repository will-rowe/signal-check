#!/usr/bin/env nextflow
pipelineVersion = '0.1'

/*
    check some stuff
*/
if (params.reads == '') {
    exit 1, "Please specify the read data (use --reads)"
}
if (params.assembler == '') {
        exit 1, "Please specify an assembler (use --assembler)'"
}
if (params.assembler != 'miniasm' && params.assembler != 'redbean') {
        exit 1, "The assembler must be 'miniasm' or 'redbean'. You tried: --assembler ${params.assembler}'"
}
if (params.output == '') {
    exit 1, "Please specify an output directory for the results (use --output)"
}

/*
    print some stuff
*/
log.info "========================================="
log.info " pipeline v${pipelineVersion}"
log.info "========================================="

def summary = [:]
summary['Reads']        = params.reads
summary['Fast5 dir']        = params.fast5
summary['Racon iterat.'] = params.raconRounds
summary['Medaka model'] = params.medakaModel
summary['Output dir']   = params.output
summary['Working dir']  = workflow.workDir
summary['Max. memory']   = params.mem
summary['Max. CPUs']     = params.cpus
summary['Profile'] = workflow.profile
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Completition message
workflow.onComplete {
    log.info "========================================="
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

/*
    do some trimmingAdapters with porechop
*/
process trimmingAdapters {
    input:
	   file(reads) from file(params.reads)

    output:
	   file('reads.trimmed.fastq') into trimmed_reads

	script:
        """
    	porechop -i "${reads}" -t "${task.cpus}" -o reads.trimmed.fastq
        """
}

// duplicate the trimmed_read channel so that it can be used by multiple downstream processes
trimmed_reads.into { trimmed_reads_for_assembly; trimmed_reads_for_polishing }

/*
    do an assembly with miniasm or redbean
*/
process assemblingReads {
    publishDir params.output, mode: 'copy', pattern: 'assembly-unpolished.fasta'
    echo true

    input:
        file reads from trimmed_reads_for_assembly

    output:
        file('assembly-unpolished.fasta') into assembly

    script:
        if(params.assembler == 'miniasm')
            """
            echo "[info] using miniasm for assembly"
            minimap2 -x ava-ont -F 200 -t "${task.cpus}" "${reads}" "${reads}" > "${reads}.paf"
            miniasm -e2 -n1 -s 500 -R -f "${reads}" "${reads}.paf" > "${reads}.gfa"
            awk '/^S/{print ">"\$2"\\n"\$3}' "${reads}.gfa" | fold > assembly-unpolished.fasta
            """
        else if(params.assembler == 'redbean')
            """
            echo "[info] using readbean for assembly (not impl. yet)"

            echo "[info] using miniasm for assembly"
            minimap2 -x ava-ont -F 200 -t "${task.cpus}" "${reads}" "${reads}" > "${reads}.paf"
            miniasm -e2 -n1 -s 500 -R -f "${reads}" "${reads}.paf" > "${reads}.gfa"
            awk '/^S/{print ">"\$2"\\n"\$3}' "${reads}.gfa" | fold > assembly-unpolished.fasta
            """
}

// copy the assembly channel so that it can be used by multiple downstream processes
assembly.into { assembly_for_polishing; assembly_for_signal_based_polishing; assembly_for_read_subsampling }

/*
    do some subsampling of the reads
*/
process subsamplingReads {
    input:
        file(assembly) from assembly_for_read_subsampling
        file(reads) from file(params.reads)

    output:
        file('sub_sampled.bam') into subsampled_bam
        file('sub_sampled.reads.fq') into subsampled_reads

    script:
        """
        minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o reads.sorted.bam -T reads.tmp -
        samtools index reads.sorted.bam
        subsample_bam --proportional -o sub_sampled reads.sorted.bam "${params.subSamplingDepth}"
        mv sub_sampled*.bam sub_sampled.bam
        samtools fastq sub_sampled.bam > sub_sampled.reads.fq
        """
}

/*
    do some polishing without using the signal (X rounds racon, Y round medaka)
*/
process polishingWithoutSignal {
	publishDir params.output, mode: 'copy', pattern: 'assembly-polished-without-using-signal.fasta'
    echo true
    input:
        file(assembly) from assembly_for_polishing
    	file(reads) from trimmed_reads_for_polishing

    output:
	   file('assembly-polished-without-using-signal.fasta') into assembly_polished_without_using_signal

	script:
        """
        minimap2 -x map-ont -t "${task.cpus}" "${assembly}" "${reads}" > assembly-racon1.paf
        racon -t "${task.cpus}" "${reads}" assembly-racon1.paf "${assembly}" > assembly-racon1.fasta

        remainingRounds="\$((${params.raconRounds}-1))"
        for (( i=1; i<=\$remainingRounds; i++ ))
        do
            ii=\$(( \$i + 1 ))
            minimap2 -x map-ont -t "${task.cpus}" assembly-racon\$i.fasta "${reads}" > assembly-racon\$ii.paf
            racon -t "${task.cpus}" "${reads}" assembly-racon\$ii.paf assembly-racon\$i.fasta > assembly-racon\$ii.fasta
        done

        medaka_consensus -i "${reads}" -d assembly-racon4.fasta -o medaka -m "${params.medakaModel}" -t "${task.cpus}"

        cp medaka/consensus.fasta assembly-polished-without-using-signal.fasta
        """
}

/*
    do some polishing using the signal
*/
process polishingWithSignal {
    publishDir params.output, mode: 'copy', pattern: 'assembly-polished-using-signal.fasta'
    if (params.fast5 == '')
        echo true

    input:
        file(assembly) from assembly_for_signal_based_polishing
        file(reads) from subsampled_reads
        file(bam) from subsampled_bam
        val(fast5_dir) from params.fast5

    output:
        file('assembly-polished-using-signal.fasta') into assembly_polished_using_signal

    script:
        if (params.fast5 != '')
            """
            nanopolish index -d "${fast5_dir}" "${reads}"
            samtools index "${bam}"
            nanopolish variants --consensus -t "${task.cpus}" -r "${reads}" -b "${bam}" -g "${assembly}" -o polished.vcf
            nanopolish vcf2fasta --skip-checks -g "${assembly}" polished.vcf > assembly-polished-using-signal.fasta;


            #nanopolish_makerange.py "${assembly}" | parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus -t 1 -w {1} -r "${reads}" -b reads.sorted.bam -g "${assembly}" -o polished.{1}.vcf
            #for i in polished.*.vcf
            #do;
            #base=\${i%%.vcf};
            #nanopolish vcf2fasta --skip-checks -g "${assembly}" \$i > \${base}.fasta;
            #done;
            #nanopolish_merge.py polished.*.fasta > assembly-polished-using-signal.fasta
            """
        else
            """
            echo "[warn] can't do signal-based polishing without signal"
            echo "[warn] skipping signal-based polishing and touching dummy outfile"
            touch assembly-polished-using-signal.fasta
            """
}

/*
    do some assessment of assemblies
*/
process assessAssemblies {
    publishDir params.output, mode: 'copy', pattern: 'quast_reports.tar'

    input:
        file(assembly1) from assembly_polished_without_using_signal
        file(assembly2) from assembly_polished_using_signal

    output:
        file('quast_reports.tar') into assembly_assessed

    script:
        """
        quast.py -o quast_polished_without_signal -t "${task.cpus}" --circos "${assembly1}"
        quast.py -o quast_polished_with_signal -t "${task.cpus}" --circos "${assembly2}"

        mkdir quast_reports
        mv quast_polished_without_signal quast_reports/
        mv quast_polished_with_signal quast_reports/
        tar -cvf quast_reports.tar quast_reports
        """
}

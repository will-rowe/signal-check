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

/*
    do some adapterTrimming with porechop
*/
process adapterTrimming {
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
process assembly {
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
            minimap2 -x ava-ont -t "${task.cpus}" "${reads}" "${reads}" > "${reads}.paf"
            miniasm -f "${reads}" "${reads}.paf" > "${reads}.gfa"
            awk '/^S/{print ">"\$2"\\n"\$3}' "${reads}.gfa" | fold > assembly-unpolished.fasta
            """
        else if(params.assembler == 'redbean')
            """
            echo "[info] using readbean for assembly (not impl. yet)"


            minimap2 -x ava-ont -t "${task.cpus}" "${reads}" "${reads}" > "${reads}.paf"
            miniasm -f "${reads}" "${reads}.paf" > "${reads}.gfa"
            awk '/^S/{print ">"\$2"\\n"\$3}' "${reads}.gfa" | fold > assembly-unpolished.fasta
            """
}

// duplicate the assembly channel so that it can be used by multiple downstream processes
assembly.into { assembly_for_polishing; assembly_for_signal_based_polishing }

/*
    do some polishing without using the signal (4 rounds racon, 1 round medaka)
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
        file(reads) from file(params.reads)
        val(fast5_dir) from params.fast5

    output:
        file('assembly-polished-using-signal.fasta') into assembly_polished_using_signal

    script:
        if (params.fast5 != '')
            """
            nanopolish index -d "${fast5_dir}" "${reads}"

            minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o reads.sorted.bam -T reads.tmp -
            samtools index reads.sorted.bam
            nanopolish_makerange.py "${assembly}" | \
            parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus polished.{1}.fa -w {1} -r "${reads}" -b reads.sorted.bam -g "${assembly}" -t 1 --min-candidate-frequency 0.1
            nanopolish_merge.py polished.*.fa > assembly-polished-using-signal.fasta
            """
        else
            """
            echo "[warn] can't do signal-based polishing without signal"
            echo "[warn] skipping signal-based polishing and touching dummy outfile"
            touch assembly-polished-using-signal.fasta
            """
}

#!/usr/bin/env nextflow
pipelineVersion = '0.1'

/*
    check some stuff
*/
if (params.fastqDir == '') {
    exit 1, "Please specify the directory containing the fastq data - this is usually the fastq_pass directory in the MinKNOW output (use --fastqDir)"
}
if (params.fast5Dir == '') {
    exit 1, "Please specify the directory containing the fast5 data - this is usually the fast5_pass directory in the MinKNOW output (use --fast5Dir)"
}
if (params.barcodes.size() == 0) {
    exit 1, "Please specify the barcodes to use (e.g --barcodes 09,10,11)"
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
log.info "-------------------------------------------------------"
log.info " de novo long read assembly pipeline v${pipelineVersion}"
log.info "-------------------------------------------------------"

def summary = [:]
summary['fastq dir']        = params.fastqDir
summary['fast5 dir']        = params.fast5Dir
summary['Barcodes'] = params.barcodes
summary['Sequencing kit'] = params.seqKit
if (params.subSamplingDepth != 0) {
    summary['Sampling depth'] = params.subSamplingDepth
} else {
    summary['Sampling depth'] = "na"
}
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
log.info "-------------------------------------------------------"

// completion message
workflow.onComplete {
    log.info "-------------------------------------------------------"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

/*
    do some demuxing and trimming
*/
process demuxingReads {
    publishDir params.output, mode: 'copy', pattern: 'demuxed_reads/barcode-*.fastq'
    echo false

    input:
       val(fastqDir) from params.fastqDir
       

    output:
	   file('demuxed_reads/barcode-*.fastq') into trimmed_reads_for_assembly

	script:
        """
        echo "[info] demuxing reads, trimming and keeping barcodes "${params.barcodes}""
        cat "${fastqDir}"/*.fastq | qcat --threads "${task.cpus}" --guppy --trim --kit "${params.seqKit}" --min-score "${params.minQualScore}" -b demuxed_reads
        qcat-parser.py demuxed_reads "${params.barcodes}"
        """
}

/*
    do some assemblies with miniasm or redbean
*/
process assemblingReads {
    publishDir params.output, mode: 'copy', pattern: '*.assembly-unpolished.fasta'
    echo false

    input:
        file(reads) from trimmed_reads_for_assembly.flatten()

    output:
        file('*.assembly-unpolished.fasta') into assembly_for_correction
        file(reads) into trimmed_reads_for_correction

    script:
        if(params.assembler == 'miniasm')
            """
            echo "[info] using miniasm for assembly of "${reads}""
            minimap2 -x ava-ont -F 200 -t "${task.cpus}" "${reads}" "${reads}" > "${reads.getBaseName()}.paf"
            miniasm -e2 -n1 -s 500 -R -f "${reads}" "${reads.getBaseName()}.paf" > "${reads.getBaseName()}.gfa"
            awk '/^S/{print ">"\$2"\\n"\$3}' "${reads.getBaseName()}.gfa" | fold > ${reads.getBaseName()}.assembly-unpolished.fasta
            """
        else if(params.assembler == 'redbean')
            """
            echo "[info] using readbean for assembly of "${reads}""
            wtdbg2 -x ont -g 20Kb -i "${reads}" -t "${task.cpus}" -fo "${reads.getBaseName()}"
            wtpoa-cns -t "${task.cpus}" -i "${reads.getBaseName()}".ctg.lay.gz -fo "${reads.getBaseName()}".assembly-unpolished.fasta
            """
}

/*
    do some assembly correction with x rounds of racon
*/
process correctingAssemblyWithRacon {
    input:
        file(assembly) from assembly_for_correction
    	file(reads) from trimmed_reads_for_correction

    output:
	   file('*.assembly.racon-corrected.fasta') into corrected_assembly
       file(reads) into trimmed_reads

	script:
        """
        minimap2 -x map-ont -t "${task.cpus}" "${assembly}" "${reads}" > ${reads.getBaseName()}-racon1.paf
        racon -t "${task.cpus}" "${reads}" ${reads.getBaseName()}-racon1.paf "${assembly}" > ${reads.getBaseName()}-racon1.fasta

        remainingRounds="\$((${params.raconRounds}-1))"
        for (( i=1; i<=\$remainingRounds; i++ ))
        do
            ii=\$(( \$i + 1 ))
            minimap2 -x map-ont -t "${task.cpus}" ${reads.getBaseName()}-racon\$i.fasta "${reads}" > ${reads.getBaseName()}-racon\$ii.paf
            racon -t "${task.cpus}" "${reads}" ${reads.getBaseName()}-racon\$ii.paf ${reads.getBaseName()}-racon\$i.fasta > ${reads.getBaseName()}-racon\$ii.fasta
            rm ${reads.getBaseName()}-racon\$i.fasta
        done

        mv ${reads.getBaseName()}-racon*.fasta ${reads.getBaseName()}.assembly.racon-corrected.fasta
        """
}


/*
    do some read subsampling if requested
*/
process subsamplingReads {
    echo false
    
    input:
    	file(assembly) from corrected_assembly
        file(reads) from trimmed_reads

    output:
        file(assembly) into subsampled_assembly
        file(reads) into subsampled_reads_original
        file('*.sub-sampled.fq') into subsampled_reads

    script:
        if (params.subSamplingDepth != 0)
            """
            echo "[info] subsampling reads for ${reads.getBaseName()}, using "${reads}" and "${assembly}""
            minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o ${reads.getBaseName()}.assembly-alignment.bam -
            samtools index ${reads.getBaseName()}.assembly-alignment.bam
            subsample_bam -t "${task.cpus}" -o sub_sampled ${reads.getBaseName()}.assembly-alignment.bam "${params.subSamplingDepth}"
            mv sub_sampled*.bam ${reads.getBaseName()}.assembly-alignment.bam
            samtools fastq ${reads.getBaseName()}.assembly-alignment.bam > ${reads.getBaseName()}.sub-sampled.fq
            """
        else
            """
            echo "[info] skipping read sub-sampling for ${reads.getBaseName()}"
            mv "${reads}" ${reads.getBaseName()}.sub-sampled.fq
            """
}


// copy the assembly and reads channel so that we can fork the pipeline
subsampled_assembly.into {assembly_for_medaka_f1; assembly_for_nanopolish_f2}
subsampled_reads_original.into {reads_for_medaka_f1; reads_for_medaka_f2}
subsampled_reads.into {reads_for_nanopolish_f1; reads_for_nanopolish_f2}

/*
    FORK 1:
*/

/*
    do some polishing with medaka
*/
process polishingWithMedaka {
	publishDir params.output, mode: 'copy', pattern: '*.assembly-corrected.medaka-polished.fasta'

    input:
        file(assembly) from assembly_for_medaka_f1
    	file(reads) from reads_for_medaka_f1
        file(subsampledReads) from reads_for_nanopolish_f1

    output:
	   file('*.assembly-corrected.medaka-polished.fasta') into medaka_polished_assembly
       file(subsampledReads) into medaka_reads_for_repolishing

	script:
        """
        medaka_consensus -i "${reads}" -d "${assembly}" -o medaka -m "${params.medakaModel}" -t "${task.cpus}"
        awk '/^>/{print ">contig" ++i; next}{print}' < medaka/consensus.fasta > ${reads.getBaseName()}.assembly-corrected.medaka-polished.fasta
        """
}

medaka_polished_assembly.into {assembly_polished_without_using_signal; medaka_polished_assembly_for_repolishing}

/*
    do some further polishing of the medaka polished assembly, using nanopolish
*/
process repolishingWithNanopolish {
	publishDir params.output, mode: 'copy', pattern: '*.assembly-corrected.medaka-polished.nanopolish-repolished.fasta'

    input:
        file(assembly) from medaka_polished_assembly_for_repolishing
    	file(reads) from medaka_reads_for_repolishing
        val(fast5Dir) from params.fast5Dir

    output:
	   file('*.assembly-corrected.medaka-polished.nanopolish-repolished.fasta') into assembly_medaka_then_nanopolish

    script:
        """
        minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o ${reads.getBaseName()}.assembly-alignment.bam -
        samtools index ${reads.getBaseName()}.assembly-alignment.bam

        nanopolish index -d "${fast5Dir}" "${reads}"
        nanopolish variants --consensus -t "${task.cpus}" -r "${reads}" -b ${reads.getBaseName()}.assembly-alignment.bam -g "${assembly}" -o polished.vcf
        nanopolish vcf2fasta --skip-checks -g "${assembly}" polished.vcf > "${reads.getBaseName()}".assembly-corrected.medaka-polished.nanopolish-repolished.fasta;
        """
}

/*
    FORK 2:
*/

/*
    do some polishing with nanopolish
*/
process polishingWithNanopolish {
    publishDir params.output, mode: 'copy', pattern: '*.assembly-corrected.nanopolish-polished.fasta'

    input:
        file(assembly) from assembly_for_nanopolish_f2
        file(reads) from reads_for_nanopolish_f2
        file(origReads) from reads_for_medaka_f2
        val(fast5Dir) from params.fast5Dir

    output:
        file('*.assembly-corrected.nanopolish-polished.fasta') into nanopolished_assembly
        file(origReads) into reads_for_repolishing

    script:
        """
        minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o ${reads.getBaseName()}.assembly-alignment.bam -
        samtools index ${reads.getBaseName()}.assembly-alignment.bam

        nanopolish index -d "${fast5Dir}" "${reads}"
        nanopolish variants --consensus -t "${task.cpus}" -r "${reads}" -b "${reads.getBaseName()}.assembly-alignment.bam" -g "${assembly}" -o polished.vcf
        nanopolish vcf2fasta --skip-checks -g "${assembly}" polished.vcf > "${reads.getBaseName()}".assembly-corrected.nanopolish-polished.fasta;
        """
}

nanopolished_assembly.into {assembly_polished_using_signal; nanopolished_assembly_for_repolishing}

/*
    do some further polishing of the nanopolished assembly, using medaka and the original reads
*/
process repolishingWithMedaka {
	publishDir params.output, mode: 'copy', pattern: '*.assembly-corrected.nanopolish-polished.medaka-repolished.fasta'

    input:
        file(assembly) from nanopolished_assembly_for_repolishing
    	file(reads) from reads_for_repolishing

    output:
	   file('*.assembly-corrected.nanopolish-polished.medaka-repolished.fasta') into assembly_nanopolish_then_medaka

	script:
        """
        medaka_consensus -i "${reads}" -d "${assembly}" -o medaka -m "${params.medakaModel}" -t "${task.cpus}"
        awk '/^>/{print ">contig" ++i; next}{print}' < medaka/consensus.fasta > ${reads.getBaseName()}.assembly-corrected.nanopolish-polished.medaka-repolished.fasta
        """
}

/*
    REJOIN FORKS:
*/

/*
    do some assessment of assemblies
*/
process assessAssemblies {
    publishDir params.output, mode: 'copy', pattern: 'quast_reports.tar'

    input:
        file(assembly1) from assembly_polished_without_using_signal
        file(assembly2) from assembly_polished_using_signal
        file(assembly3) from assembly_nanopolish_then_medaka
        file(assembly4) from assembly_medaka_then_nanopolish

    output:
        file('quast_reports.tar') into assembly_assessed

    script:
        """
        quast.py -o quast_polished_with_medaka -t "${task.cpus}" --circos "${assembly1}"
        quast.py -o quast_polished_with_nanopolish -t "${task.cpus}" --circos "${assembly2}"
        quast.py -o quast_nanopolish_then_medaka -t "${task.cpus}" --circos "${assembly3}"
        quast.py -o quast_medaka_then_nanopolish -t "${task.cpus}" --circos "${assembly4}"

        mkdir quast_reports
        mv quast_polished_with_medaka quast_reports/
        mv quast_polished_with_nanopolish quast_reports/
        mv quast_nanopolish_then_medaka quast_reports/
        mv quast_medaka_then_nanopolish quast_reports/
        tar -cvf quast_reports.tar quast_reports
        """
}

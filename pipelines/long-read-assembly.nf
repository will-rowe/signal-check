#!/usr/bin/env nextflow
pipelineVersion = '0.1'

/*
    check some stuff
*/
if (params.inputDir == '') {
    exit 1, "Please specify the input directory - this is the MinKNOW output dir containing fastq_pass and fast5_pass (use --inputDir)"
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
log.info " long read assembly pipeline v${pipelineVersion}"
log.info "-------------------------------------------------------"

def summary = [:]
summary['Input dir']        = params.inputDir
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
    echo true

    input:
       val(inputDir) from params.inputDir
       

    output:
	   file('demuxed_reads/barcode-*.fastq') into trimmed_reads_for_assembly

	script:
        """
        echo "[info] demuxing reads, trimming and keeping barcodes "${params.barcodes}""
        cat "${inputDir}"/fastq_pass/*.fastq | qcat --threads "${task.cpus}" --guppy --trim --kit "${params.seqKit}" --min-score "${params.minQualScore}" -b demuxed_reads
        qcat-parser.py demuxed_reads "${params.barcodes}"
        """
}

/*
    do some assemblies with miniasm or redbean
*/
process assemblingReads {
    publishDir params.output, mode: 'copy', pattern: '*.assembly-unpolished.fasta'
    echo true

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
            wtdbg2 -x ont -i "${reads}" -t "${task.cpus}" -fo "${reads.getBaseName()}"
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
	   file('*.assembly-corrected.fasta') into corrected_assembly
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

        mv ${reads.getBaseName()}-racon*.fasta ${reads.getBaseName()}.assembly-corrected.fasta
        """
}

// copy the assembly and reads channel so that we can fork the pipeline
corrected_assembly.into {assembly_for_subsampling; assembly_for_medaka}
trimmed_reads.into {reads_for_subsampling; reads_for_medaka}

/*
    do some polishing with medaka
*/
process polishingWithMedaka {
	publishDir params.output, mode: 'copy', pattern: '*.assembly-corrected.medaka-polished.fasta'
    echo false

    input:
        file(assembly) from assembly_for_medaka
    	file(reads) from reads_for_medaka

    output:
	   file('*.assembly-corrected.medaka-polished.fasta') into assembly_polished_without_using_signal

	script:
        """
        echo "[info] medaka for correction of "${assembly}" using "${reads}""
        medaka_consensus -i "${reads}" -d "${assembly}" -o medaka -m "${params.medakaModel}" -t "${task.cpus}"
        cp medaka/consensus.fasta ${reads.getBaseName()}.assembly-corrected.medaka-polished.fasta
        """
}

/*
    do some optional subsampling of the reads prior to nanopolish
*/
process subsamplingReads {
    echo true

    input:
        file(assembly) from assembly_for_subsampling
        file(reads) from reads_for_subsampling

    output:
        file(assembly) into assembly_for_nanopolish
        file('*.sub-sampled.fq') into reads_for_nanopolish
        file('*.assembly-alignment.bam') into alignment_for_nanopolish

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
            minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o ${reads.getBaseName()}.assembly-alignment.bam -
            samtools index ${reads.getBaseName()}.assembly-alignment.bam
            mv "${reads}" ${reads.getBaseName()}.sub-sampled.fq
            """
}

/*
    do some polishing with nanopolish
*/
process polishingWithNanopolish {
    publishDir params.output, mode: 'copy', pattern: '*.assembly-corrected.nanopolish-polished.fasta'

    input:
        file(assembly) from assembly_for_nanopolish
        file(reads) from reads_for_nanopolish
        file(bam) from alignment_for_nanopolish
        val(inputDir) from params.inputDir

    output:
        file('*.assembly-corrected.nanopolish-polished.fasta') into assembly_polished_using_signal

    script:
        """
        nanopolish index -d "${inputDir}" "${reads}"
        samtools index "${bam}"
        nanopolish variants --consensus -t "${task.cpus}" -r "${reads}" -b "${bam}" -g "${assembly}" -o polished.vcf
        nanopolish vcf2fasta --skip-checks -g "${assembly}" polished.vcf > "${reads.getBaseName()}".assembly-corrected.nanopolish-polished.fasta;

        #nanopolish_makerange.py "${assembly}" | parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus -t 1 -w {1} -r "${reads}" -b "${bam}" -g "${assembly}" -o polished.{1}.vcf
        #do;
        #for i in polished.*.vcf
        #base=\${i%%.vcf};
        #nanopolish vcf2fasta --skip-checks -g "${assembly}" \$i > \${base}.fasta;
        #done;
        #nanopolish_merge.py polished.*.fasta > "${reads.getBaseName()}".assembly-corrected.nanopolish-polished.fasta;
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

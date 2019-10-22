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
if (params.refGenomes == '') {
    exit 1, "Please specify the reference genomes multifasta (use --refGenomes)" 
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
log.info "long read assembly pipeline v${pipelineVersion}"
log.info "-------------------------------------------------------"

def summary = [:]
summary['Directory with fastq']        = params.fastqDir
summary['Directory with fast5']        = params.fast5Dir
summary['Reference genomes'] = params.refGenomes
if (params.label != '') {
    summary['File prepend label'] = params.label
}
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

// error message
workflow.onError {
    log.info "-------------------------------------------------------"
    log.info "Oops... Pipeline execution stopped with the following message:"
    println "${workflow.errorMessage}"
}

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    START THE PROCESSES:
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    do some demuxing and trimming
*/
process demuxingReads {
    publishDir params.output, mode: 'copy', pattern: 'demuxed_reads/*arcode*.fastq'
    echo false

    input:
       val(fastqDir) from params.fastqDir

    output:
	   file('demuxed_reads/*arcode-*.fastq') into trimmed_reads

	script:
        """
        echo "[info] demuxing reads, trimming and keeping barcodes "${params.barcodes}""
        cat "${fastqDir}"/*.fastq | qcat --threads "${task.cpus}" --guppy --trim --kit "${params.seqKit}" --min-score "${params.minQualScore}" -b demuxed_reads
        qcat-parser.py demuxed_reads "${params.barcodes}" "${params.label}"
        """
}

trimmed_reads.into {reads_for_alignment; reads_for_rg_assembly; reads_for_dn_assembly; reads_for_quast}

/*
    do an alignment to the specified reference genome
*/
process referenceAlignment {
    publishDir params.output + "/reference-alignment", mode: 'copy', pattern: '*.png'
    errorStrategy 'terminate'
    echo false

    input:
        file(reads) from reads_for_alignment.flatten()
        file(refGenomes) from params.refGenomes

    output:
        file('*.ref-alignment.bam') into reference_alignments
        file('refGenome.fasta') into reference_sequence
        file(reads) into reads_for_subsampling
        file('*.png') into alignment_pngs

    script:
        """
        get-ref.py ${reads} ${params.refGenomes} refGenome.fasta
        minimap2 -ax map-ont -t "${task.cpus}" refGenome.fasta "${reads}" | \
        samtools view -bS - | \
        samtools sort - -o ${reads.getBaseName()}.ref-alignment.bam
        samtools index ${reads.getBaseName()}.ref-alignment.bam

        # check that reads map back to assembly
        mappedReads=\$(samtools view -F 4 -c ${reads.getBaseName()}.ref-alignment.bam)
        echo \$mappedReads
        if [ \$mappedReads == 0 ]; then
            echo "no reads mapped to the reference - wrong reference?"
            exit 1
        fi

        samtools depth -a ${reads.getBaseName()}.ref-alignment.bam > ${reads.getBaseName()}.ref-alignment.depth        
        plot-bam-depth.py ${reads.getBaseName()}.ref-alignment.depth ${reads.getBaseName()}.ref-alignment
        """
}

reference_alignments.into{alignments_for_subsampling; alignments_for_variant_calling}
reference_sequence.into{reference_for_medaka_variant_calling; reference_for_nanopolish_variant_calling}

/*
    do some read subsampling if requested
*/
process subsamplingReads {
    echo false
    
    input:
    	file(alignment) from alignments_for_subsampling
        file(reads) from reads_for_subsampling

    output:
        file('*.sub-sampled.fq') into subsampled_reads

    script:
        if (params.subSamplingDepth != 0)
            """
            samtools index ${alignment}
            subsample_bam -t "${task.cpus}" -o sub_sampled "${alignment}" "${params.subSamplingDepth}"
            for i in sub_sampled*.bam; do
            samtools fastq \$i >> ${reads.getBaseName()}.sub-sampled.fq;
            done
            """
        else
            """
            echo "[info] skipping read sub-sampling for ${reads.getBaseName()}"
            mv "${reads}" ${reads.getBaseName()}.sub-sampled.fq
            """
}

/*
    run nanopolish index on the subsampled reads (or original if no subsampling was requested)
*/
process nanopolishIndexing {
    input:
        file(reads) from subsampled_reads
        val(fast5Dir) from params.fast5Dir

    output:
        set file(reads), file('*.index'), file('*.fai'), file('*.gzi'), file('*.readdb') into nanopolish_index_files

    script:
        """
        nanopolish index -d "${fast5Dir}" "${reads}"
        """
}

// copy the nanopolish index files for each fork
nanopolish_index_files.into {nanopolish_index_files_for_dn_assemblies; nanopolish_index_files_for_variant_calling}

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    SPLIT THE PIPELINE: variant calling, reference guided assembly, de novo assembly
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    VARIANT CALLING:
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    variant call with medaka
*/
process variantcallWithMedaka {
    publishDir params.output + "/reference-alignment", mode: 'copy', pattern: '*.medaka-consensus.fasta'

    input:
        file(alignment) from alignments_for_variant_calling
        file(reference) from reference_for_medaka_variant_calling

    output:
        file('*.medaka-consensus.fasta') into medaka_consensus

    script:
        """
        medaka_variant -f ${reference} -i ${alignment} -m ${params.medakaModel} -s ${params.medakaModel} -d -t ${task.cpus} -o ${alignment.getBaseName()}
        mv ${alignment.getBaseName()}/round_1_phased.vcf ${alignment.getBaseName()}.vcf
        bgzip ${alignment.getBaseName()}.vcf
        tabix ${alignment.getBaseName()}.vcf.gz
        bcftools consensus -i '%QUAL>=20' -f ${reference} ${alignment.getBaseName()}.vcf.gz > ${alignment.getBaseName()}.medaka-consensus.fasta
        """
}

/*
    variant call with nanopolish
*/
process variantcallWithNanopolish {
    publishDir params.output + "/reference-alignment", mode: 'copy', pattern: '*.nanopolish-consensus.fasta'

    input:
        file(reference) from reference_for_nanopolish_variant_calling
        set file(reads), file(index), file(fai), file(gzi), file(readdb) from nanopolish_index_files_for_variant_calling

    output:
        file('*.nanopolish-consensus.fasta') into nanopolish_consensus

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} | \
        samtools view -bS - | \
        samtools sort - -o ${reads.getBaseName()}.ref-alignment.bam
        samtools index ${reads.getBaseName()}.ref-alignment.bam
        nanopolish variants --reads ${reads} --bam ${reads.getBaseName()}.ref-alignment.bam --genome ${reference} -t ${task.cpus} --ploidy 1 --snps -o nanopolish_variants.vcf
        bgzip nanopolish_variants.vcf
        tabix nanopolish_variants.vcf.gz
        bcftools consensus -i '%QUAL>=20' -f ${reference} nanopolish_variants.vcf.gz > ${reads.getBaseName()}.nanopolish-consensus.fasta

        #margin_cons.py ${reference} nanopolish_variants.vcf ${reads.getBaseName()}.ref-alignment.bam > ${reads.getBaseName()}.nanopolish-consensus.fasta
        """
}

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    REFERENCE GUIDED ASSEMBLY:
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    kick of a reference guided assembly
*/
process assemblingReadsRG {
    publishDir params.output + "/reference-guided-assembly", mode: 'copy', pattern: '*.rg-assembly.racon*'
    errorStrategy 'terminate'
    echo false

    input:
        file(reads) from reads_for_rg_assembly.flatten()
        file(refGenomes) from params.refGenomes

    output:
        file('*.rg-assembly.racon.fasta') into rg_assemblies
        file('*.png') into rg_pngs

    script:
        """
        get-ref.py ${reads} ${params.refGenomes} refGenome.fasta

        mini_assemble -i "${reads}" -r refGenome.fasta -m "${params.raconRounds}" -t "${task.cpus}" -o pomoxis -p ${reads.getBaseName()}.assembly.racon
        mv pomoxis/${reads.getBaseName()}.assembly.racon_final.fa ${reads.getBaseName()}.rg-assembly.racon.fasta
        minimap2 -ax map-ont -t "${task.cpus}" ${reads.getBaseName()}.rg-assembly.racon.fasta "${reads}" | \
        samtools view -bS - | \
        samtools sort - -o ${reads.getBaseName()}.rg-assembly.racon.bam
        samtools index ${reads.getBaseName()}.rg-assembly.racon.bam

        # check that reads map back to assembly
        mappedReads=\$(samtools view -F 4 -c ${reads.getBaseName()}.rg-assembly.racon.bam)
        echo \$mappedReads
        if [ \$mappedReads == 0 ]; then
            echo "no reads mapped back to the reference guided assembly - wrong reference?"
            exit 1
        fi

        samtools depth -a ${reads.getBaseName()}.rg-assembly.racon.bam > ${reads.getBaseName()}.rg-assembly.racon.depth        
        plot-bam-depth.py ${reads.getBaseName()}.rg-assembly.racon.depth ${reads.getBaseName()}.rg-assembly.racon
        """
}

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    DE NOVO ASSEMBLY:
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    do some assemblies with miniasm or redbean
*/
process assemblingReadsDN {
    publishDir params.output + "/de-novo-assembly", mode: 'copy', pattern: '*.dn-assembly.raw.fasta'
    echo false

    input:
        file(reads) from reads_for_dn_assembly.flatten()

    output:
        file('*.dn-assembly.raw.fasta') into dn_assembly_for_correction
        file(reads) into trimmed_reads_for_correction

    script:
        if(params.assembler == 'miniasm')
            """
            echo "[info] using miniasm for assembly of "${reads}""
            minimap2 -x ava-ont -F 200 -t "${task.cpus}" "${reads}" "${reads}" > "${reads.getBaseName()}.paf"
            miniasm -e2 -n1 -s 500 -R -f "${reads}" "${reads.getBaseName()}.paf" > "${reads.getBaseName()}.gfa"
            awk '/^S/{print ">"\$2"\\n"\$3}' "${reads.getBaseName()}.gfa" | fold > assembly.fasta
            awk '/^>/{print ">assembly.raw.contig" ++i; next}{print}' < assembly.fasta > "${reads.getBaseName()}".dn-assembly.raw.fasta
            """
        else if(params.assembler == 'redbean')
            """
            echo "[info] using readbean for assembly of "${reads}""
            wtdbg2 -x ont -g 20Kb -i "${reads}" -t "${task.cpus}" -fo "${reads.getBaseName()}"
            wtpoa-cns -t "${task.cpus}" -i "${reads.getBaseName()}".ctg.lay.gz -fo assembly.fasta
            awk '/^>/{print ">assembly.raw.contig" ++i; next}{print}' < assembly.fasta > "${reads.getBaseName()}".dn-assembly.raw.fasta
            """
}

/*
    do some assembly correction with x rounds of racon
*/
process correctingAssemblyWithRaconDN {
    input:
        file(assembly) from dn_assembly_for_correction
    	file(reads) from trimmed_reads_for_correction

    output:
	   file('*.dn-assembly.racon.fasta') into corrected_dn_assembly
       file(reads) into trimmed_reads_post_racon

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

        mv ${reads.getBaseName()}-racon*.fasta ${reads.getBaseName()}.dn-assembly.racon.fasta
        """
}

/*
    fork the de novo assembly pathway
*/
corrected_dn_assembly.into{dn_assembly_for_medaka; dn_assembly_for_nanopolish}
trimmed_reads_post_racon.into{reads_for_medaka_f1; reads_for_medaka_f2}

/*
    do some polishing with medaka
*/
process polishingWithMedakaDN {
	publishDir params.output + "/de-novo-assembly", mode: 'copy', pattern: '*.dn-assembly.racon.medaka.fasta'

    input:
        file(assembly) from dn_assembly_for_medaka
    	file(reads) from reads_for_medaka_f1

    output:
	   file('*.dn-assembly.racon.medaka.fasta') into medaka_polished_dn_assembly

	script:
        """
        medaka_consensus -i "${reads}" -d "${assembly}" -o medaka -m "${params.medakaModel}" -t "${task.cpus}"
        awk '/^>/{print ">medaka.contig" ++i; next}{print}' < medaka/consensus.fasta > ${assembly.getSimpleName()}.dn-assembly.racon.medaka.fasta
        """
}

medaka_polished_dn_assembly.into{dn_assembly_polished_without_using_signal; medaka_polished_dn_assembly_for_repolishing}
nanopolish_index_files_for_dn_assemblies.into{nif_dn_repolishing; nif_dn_polishing}

/*
    do some further polishing of the medaka polished assembly, using nanopolish
*/
process repolishingWithNanopolishDN {
	publishDir params.output + "/de-novo-assembly", mode: 'copy', pattern: '*.dn-assembly.racon.medaka.nanopolish.fasta'
    input:
        file(assembly) from medaka_polished_dn_assembly_for_repolishing
        set file(reads), file(index), file(fai), file(gzi), file(readdb) from nif_dn_repolishing

    output:
	   file('*.dn-assembly.racon.medaka.nanopolish.fasta') into dn_assembly_medaka_then_nanopolish

    script:
        """
        minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o ${reads.getBaseName()}.assembly-alignment.bam -
        samtools index ${reads.getBaseName()}.assembly-alignment.bam
        nanopolish_makerange.py "${assembly}" | parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r "${reads}" -b ${reads.getBaseName()}.assembly-alignment.bam -g "${assembly}" -t 4 --min-candidate-frequency 0.1
        nanopolish vcf2fasta --skip-checks -g "${assembly}" polished.*.vcf > assembly.fasta;
        awk '/^>/{print ">medaka.nanopolish.contig" ++i; next}{print}' < assembly.fasta > "${assembly.getSimpleName()}".dn-assembly.racon.medaka.nanopolish.fasta
        """
}

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    DE NOVO FORK 2:
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    do some polishing with nanopolish
*/
process polishingWithNanopolishDN {
    publishDir params.output + "/de-novo-assembly", mode: 'copy', pattern: '*.dn-assembly.racon.nanopolish.fasta'

    input:
        file(assembly) from dn_assembly_for_nanopolish
        set file(reads), file(index), file(fai), file(gzi), file(readdb) from nif_dn_polishing

    output:
        file('*.dn-assembly.racon.nanopolish.fasta') into nanopolished_dn_assembly

    script:
        """
        minimap2 -ax map-ont -t "${task.cpus}" "${assembly}" "${reads}" | samtools sort -o ${reads.getBaseName()}.assembly-alignment.bam -
        samtools index ${reads.getBaseName()}.assembly-alignment.bam

        nanopolish_makerange.py "${assembly}" | parallel --results nanopolish.results -P "${task.cpus}" nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r "${reads}" -b ${reads.getBaseName()}.assembly-alignment.bam -g "${assembly}" -t 4 --min-candidate-frequency 0.1
        nanopolish vcf2fasta --skip-checks -g "${assembly}" polished.*.vcf > assembly.fasta;
        awk '/^>/{print ">nanopolish.contig" ++i; next}{print}' < assembly.fasta > "${assembly.getSimpleName()}".dn-assembly.racon.nanopolish.fasta
        """
}

nanopolished_dn_assembly.into {dn_assembly_polished_using_signal; nanopolished_dn_assembly_for_repolishing}

/*
    do some further polishing of the nanopolished assembly, using medaka and the original reads
*/
process repolishingWithMedakaDN {
	publishDir params.output + "/de-novo-assembly", mode: 'copy', pattern: '*.dn-assembly.racon.nanopolish.medaka.fasta'

    input:
        file(assembly) from nanopolished_dn_assembly_for_repolishing
        file(reads) from reads_for_medaka_f2

    output:
	   file('*.dn-assembly.racon.nanopolish.medaka.fasta') into dn_assembly_nanopolish_then_medaka

	script:
        """
        medaka_consensus -i "${reads}" -d "${assembly}" -o medaka -m "${params.medakaModel}" -t "${task.cpus}"
        awk '/^>/{print ">nanopolish.medaka.contig" ++i; next}{print}' < medaka/consensus.fasta > ${assembly.getSimpleName()}.dn-assembly.racon.nanopolish.medaka.fasta
        """
}

/*
    /////////////////////////////////////////////////////////////////////////////////////////////
    REJOIN DE NOVO FORKS:
    /////////////////////////////////////////////////////////////////////////////////////////////
*/

/*
    do some assessment of assemblies with quast, map reads back to with minimap2, get depth plots
    
    TODO: the assemblies can come in any order in each channel -- needs a way to pull out the reads for each assembly
*/
process assessAssemblies {
    publishDir params.output + "/de-novo-assembly", mode: 'copy', pattern: 'assembly-qc.tar'

    input:
        file(assembly1) from dn_assembly_polished_without_using_signal
        file(assembly2) from dn_assembly_polished_using_signal
        file(assembly3) from dn_assembly_nanopolish_then_medaka
        file(assembly4) from dn_assembly_medaka_then_nanopolish

    output:
        file('assembly-qc.tar') into assembly_assessed

    script:
        """
        array=("${assembly1}" "${assembly2}" "${assembly3}" "${assembly4}")

        for i in \${array[@]};
        do
        quast.py -o \${i%%.fasta}.quast -t "${task.cpus}" --circos \$i
        done

        mkdir assembly-qc
        mv *.quast assembly-qc/
        tar -cvf assembly-qc.tar assembly-qc
        """
}

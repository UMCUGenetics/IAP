package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType

class VariantCaller extends QScript {
    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One or more bam files.", shortName="I")
    var bamFiles: List[File] = Nil

    @Input(doc="Output core filename.", shortName="O", required=true)
    var outputFilename: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    @Argument(doc="Minimum phred-scaled confidence to call variants", shortName="stand_call_conf", required=true)
    var standCallConf: Int = _ //30 //default: best-practices value

    @Argument(doc="Minimum phred-scaled confidence to emit variants", shortName="stand_emit_conf", required=true)
    var standEmitConf: Int = _ //10 //default: best-practices value

    // The following arguments are all optional.
    @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
    var dbsnpFile: File = _

    @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
    var targetFile: File = _

    @Argument(doc="Amount of padding (in bp) to add to each interval", shortName="ip", required=false)
    var intervalPadding: Int = 0

    @Argument(doc="Ploidy (number of chromosomes) per sample", shortName="ploidy", required=false)
    var samplePloidy: Int = 2

    @Argument(doc="", shortName="gvcf", required=false)
    var gvcf: Boolean = false

    @Argument(doc="", shortName="sexAware", required=false)
    var sexAware: Boolean = false

    def script() {
        // Define common settings for original HC, gvcf HC and sex aware gvcf HC.
	trait HC_Arguments extends HaplotypeCaller {
	    this.reference_sequence = referenceFile

	    this.scatterCount = numScatters
	    this.memoryLimit = maxMem
	    this.num_cpu_threads_per_data_thread = numCPUThreads

	    this.stand_emit_conf = standEmitConf
	    this.stand_call_conf = standCallConf
	}

	// original HC
	if (gvcf == false && sexAware == false) {
	    val haplotypeCaller = new HaplotypeCaller with HC_Arguments

	    // All required input
	    haplotypeCaller.input_file = bamFiles
	    haplotypeCaller.out = outputFilename + ".raw_variants.vcf"

	    // Optional input
	    if (dbsnpFile != null) {
		haplotypeCaller.D = dbsnpFile
	    }
		if (targetFile != null) {
		haplotypeCaller.L :+= targetFile
		haplotypeCaller.ip = intervalPadding
	    }

	    haplotypeCaller.sample_ploidy = samplePloidy

	    //add function to queue
	    add(haplotypeCaller)

	// GVCF HC
	} else if (gvcf == true && sexAware == false) {
	    var gvcfFiles : List[File] = Nil

	    // Make gvcf per bam file
	    for (bamFile <- bamFiles) {
		val haplotypeCaller = new HaplotypeCaller with HC_Arguments

		// All required input
		haplotypeCaller.input_file :+= bamFile
		haplotypeCaller.out = swapExt(bamFile, "bam", "g.vcf.gz")

		// gVCF settings
		haplotypeCaller.emitRefConfidence = ReferenceConfidenceMode.GVCF
		haplotypeCaller.variant_index_type = GATKVCFIndexType.LINEAR
		haplotypeCaller.variant_index_parameter = 128000

		// Optional input
		if (targetFile != null) {
		    haplotypeCaller.L :+= targetFile
		    haplotypeCaller.ip = intervalPadding
		}

		haplotypeCaller.sample_ploidy = samplePloidy

		//add function to queue
		gvcfFiles :+= haplotypeCaller.out
		add(haplotypeCaller)
	    }

	    //Joint genotyping
	    val genotypeGVCFs = new GenotypeGVCFs

	    genotypeGVCFs.V = gvcfFiles
	    genotypeGVCFs.reference_sequence = referenceFile
	    genotypeGVCFs.scatterCount = numScatters
	    genotypeGVCFs.num_threads = numCPUThreads

	    genotypeGVCFs.out = outputFilename + ".raw_variants.vcf"

	    // Optional input
	    if (dbsnpFile != null) {
		genotypeGVCFs.D = dbsnpFile
	    }

	    if (targetFile != null) {
		genotypeGVCFs.L :+= targetFile
		genotypeGVCFs.ip = intervalPadding
	    }
	    // Add function to queue
	    add(genotypeGVCFs)
        }
    }
}

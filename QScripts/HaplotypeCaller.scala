package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.utils.commandline.ClassType

class VariantCaller extends QScript {
    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One or more bam files.", shortName="I", required=true)
    var bamFiles: List[File] = Nil

    @Input(doc="Output core filename.", shortName="O", required=true)
    var outputFilename: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    @Argument(doc="Minimum phred-scaled confidence to call variants", shortName="stand_call_conf", required=false)
    var standCallConf: Int = 10 //default: best-practices value

    // The following arguments are all optional.
    @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
    var dbsnpFile: File = _

    @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
    var targetFile: File = _

    @Argument(doc="Amount of padding (in bp) to add to each interval", shortName="ip", required=false)
    var intervalPadding: Int = 0

    @Argument(doc="Ploidy (number of chromosomes) per sample", shortName="ploidy", required=false)
    var samplePloidy: Int = 2

    @Argument(doc="Use gvcf analysis method", shortName="gvcf", required=false)
    var gvcf: Boolean = false
    
    @Argument(doc="Exclusive upper bounds for reference confidence GQ bands", shortName="gqb", required=false)
    @ClassType(classOf[Int])
    var GVCFGQBands: Seq[Int] = Nil

    @Argument(doc="Use sex aware calling method", shortName="sexAware", required=false)
    var sexAware: Boolean = false
    
    @Argument(doc="Define sexes used during sex aware calling", shortName="sex", required=false)
    var sexes: List[String] = Nil

    def script() {
        // Define common settings for original HC, gvcf HC.
	trait HC_Arguments extends HaplotypeCaller {
	    this.reference_sequence = referenceFile

	    this.scatterCount = numScatters
	    this.memoryLimit = maxMem
	    this.num_cpu_threads_per_data_thread = numCPUThreads

	}

	// original HC
	if (gvcf == false && sexAware == false) {
	    val haplotypeCaller = new HaplotypeCaller with HC_Arguments

	    // All required input
	    haplotypeCaller.input_file = bamFiles
	    haplotypeCaller.out = outputFilename + ".raw_variants.vcf"
	    haplotypeCaller.stand_call_conf = standCallConf

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
	}

	// GVCF HC
	else if (gvcf == true && sexAware == false) {
	    var gvcfFiles : List[File] = Nil

	    // Make gvcf per bam file
	    for (bamFile <- bamFiles) {
		val haplotypeCaller = new HaplotypeCaller with HC_Arguments

		// All required input
		haplotypeCaller.input_file :+= bamFile
		haplotypeCaller.out = swapExt(bamFile, "bam", "g.vcf.gz")

		// gVCF settings
		haplotypeCaller.emitRefConfidence = ReferenceConfidenceMode.GVCF
		haplotypeCaller.GVCFGQBands = GVCFGQBands

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
	    
	    genotypeGVCFs.stand_call_conf = standCallConf

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

        // GVCF HC Sexaware HUMAN ONLY
	else if (gvcf == true && sexAware == true) {
	    var gvcfFiles : List[File] = Nil

	    // Make gvcf per bam file
	    for ((bamFile, sex) <- (bamFiles,sexes).zipped){

		// Set common settings
		trait HC_Arguments extends HaplotypeCaller {
		    this.reference_sequence = referenceFile

		    this.scatterCount = numScatters
		    this.memoryLimit = maxMem
		    this.num_cpu_threads_per_data_thread = numCPUThreads

		    this.input_file :+= bamFile

		    this.emitRefConfidence = ReferenceConfidenceMode.GVCF
		    this.GVCFGQBands = GVCFGQBands
		}

		val haplotypeCallerAutosome = new HaplotypeCaller with HC_Arguments
		val haplotypeCallerAllosome = new HaplotypeCaller with HC_Arguments

		// HUMAN AUTOSOME
		haplotypeCallerAutosome.sample_ploidy = 2
		haplotypeCallerAutosome.intervalsString = List("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT")
		haplotypeCallerAutosome.out = swapExt(bamFile, "bam", "Autosome.g.vcf.gz")
		
		// HUMAN ALLOSOME
		haplotypeCallerAllosome.out = swapExt(bamFile, "bam", "Allosome.g.vcf.gz")
		// Sex aware
		if ( sex == "male" ){
		    haplotypeCallerAllosome.intervalsString = List("X","Y")
		    haplotypeCallerAllosome.sample_ploidy = 1
		} else if ( sex == "female" ){
		    haplotypeCallerAllosome.intervalsString = List("X")
		    haplotypeCallerAllosome.sample_ploidy = 2
		}
		
		//add functions to queue
	        add(haplotypeCallerAutosome,haplotypeCallerAllosome)

		// Combine gvcfs -> probably can use CatVariants
		val combineGVCF = new CatVariants
		combineGVCF.reference = referenceFile
		combineGVCF.V :+= haplotypeCallerAutosome.out
		combineGVCF.V :+= haplotypeCallerAllosome.out
		combineGVCF.out = swapExt(bamFile, "bam", "g.vcf.gz")
		combineGVCF.assumeSorted = true
		gvcfFiles :+= combineGVCF.out

		//add function to queue
		add(combineGVCF)
		gvcfFiles :+= combineGVCF.out
	    }

	    //Joint genotyping
	    val genotypeGVCFs = new GenotypeGVCFs
	    genotypeGVCFs.V = gvcfFiles
	    genotypeGVCFs.reference_sequence = referenceFile
	    genotypeGVCFs.scatterCount = numScatters
	    genotypeGVCFs.num_threads = numCPUThreads
	    genotypeGVCFs.out = outputFilename + ".raw_variants.vcf"
	    
	    genotypeGVCFs.stand_call_conf = standCallConf
	    
	    // Optional input
	    if (dbsnpFile != null) {
		genotypeGVCFs.D = dbsnpFile
	    }

	    // Add function to queue
	    add(genotypeGVCFs)
        }
    }
}

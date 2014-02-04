/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.phonehome._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
import org.broadinstitute.sting.queue.extensions.snpeff._
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.apache.commons.io.FilenameUtils
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion


class preProcessing extends QScript {
  // Create an alias 'qscript' to be able to access variables in the IndelRealigner.
  // 'qscript' is now the same as 'IndelRealigner.this'
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R", required=true)
  var referenceFile: File = _ // _ is scala shorthand for null
  //var referenceFile: File = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta" // testing ref genome

  @Input(doc="Bam file to genotype.", shortName="I", required=true)
  var bamFile: File = _

  @Input(doc="Output core filename.", shortName="O", required=true)
  var out: File = _

  // The following arguments are all optional.
  //@Input(doc ="Input VCF file with known indels.", shortName="known", required=false)
  //var knownIndelFile: File = _
  //var knownIndelFile: File = "/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf" //test known indels

  // This trait allows us set the variables below in one place and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    // Set the memory limit to 2 gigabytes on each command.
    this.memoryLimit = 4
    //this.num_cpu_threads_per_data_thread = 10
//  this.num_threads = 8
  }

  def script() {
    // Create functions that we may run depending on options.
    val baseRecalibrate = new BaseRecalibrator with UnifiedGenotyperArguments
    val baseRecalibrateAfter = new BaseRecalibrator with UnifiedGenotyperArguments
    val analyzeCovariates = new AnalyzeCovariates with UnifiedGenotyperArguments
    val printReads = new PrintReads with UnifiedGenotyperArguments

    // Analyze covariation in dataset
    baseRecalibrate.scatterCount = 500
    baseRecalibrate.nct = 10
    baseRecalibrate.input_file:+= qscript.bamFile
    baseRecalibrate.R = referenceFile
    baseRecalibrate.knownSites :+= "/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/1000G_phase1.indels.b37.vcf"
    baseRecalibrate.knownSites :+= "/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/dbsnp_137.b37.vcf"
    baseRecalibrate.knownSites :+= "/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf"
    baseRecalibrate.out = qscript.out + ".recal_data.table"
    
    // Analyze covariation remaining after recalibration
    baseRecalibrateAfter.scatterCount = 500
    baseRecalibrateAfter.nct = 10
    baseRecalibrateAfter.input_file:+= qscript.bamFile
    baseRecalibrateAfter.R = referenceFile
    baseRecalibrateAfter.knownSites = baseRecalibrate.knownSites
    baseRecalibrateAfter.BQSR = baseRecalibrate.out
    baseRecalibrateAfter.out = qscript.out + ".post_recal_data.table"
    
    // Generate before and after plots
    analyzeCovariates.R = referenceFile
    analyzeCovariates.before = baseRecalibrate.out
    analyzeCovariates.after = baseRecalibrateAfter.out
    analyzeCovariates.plots = qscript.out + ".bqsr.pdf"
    
    // Apply recalibration to the dataset.
    printReads.scatterCount = 500
    printReads.nct = 10
    printReads.R = referenceFile
    printReads.input_file :+= qscript.bamFile
    printReads.BQSR = baseRecalibrate.out
    printReads.out = qscript.out + ".recal.bam"

    //add(targetcreator,realigner,baseRecalibrate)
    add(baseRecalibrate,baseRecalibrateAfter,analyzeCovariates,printReads)
  }
}
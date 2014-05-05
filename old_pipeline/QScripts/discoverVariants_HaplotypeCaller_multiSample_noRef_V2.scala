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


class discoverVariants extends QScript {
  // Create an alias 'qscript' to be able to access variables in the IndelRealigner.
  // 'qscript' is now the same as 'IndelRealigner.this'
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R", required=true)
  var referenceFile: File = _ // _ is scala shorthand for null
  //var referenceFile: File = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta" // testing ref genome

   @Input(doc="One or more bam files.", shortName="I")
  var bamFiles: List[File] = Nil

  @Input(doc="Output core filename.", shortName="O", required=true)
  var out: File = _

  // The following arguments are all optional.
  @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
  var dbsnpFile: File = _

  @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
  var targetFile: File = _
  
  @Argument(doc="Memory limit", shortName="mem", required=false)
  var memLimit: Int = 4

  @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=false)
  var nct: Int = 4
  
  @Argument(doc="Number of scatters", shortName="scatters", required=false)
  var scatters: Int = 2500

  // This trait allows us set the variables below in one place and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.memoryLimit = memLimit
    this.num_cpu_threads_per_data_thread = nct
  }

  def script() {
    // Create functions that we may run depending on options.
    val haplotypeCaller = new HaplotypeCaller with UnifiedGenotyperArguments

    // Call variants in sequence data
    haplotypeCaller.scatterCount = scatters

    haplotypeCaller.input_file = bamFiles
    haplotypeCaller.D = dbsnpFile

    haplotypeCaller.R = referenceFile
    if(targetFile == Nil) {
	haplotypeCaller.L :+= targetFile
    }

    haplotypeCaller.stand_emit_conf = 10
    haplotypeCaller.stand_call_conf = 30
    haplotypeCaller.out = qscript.out + ".raw_variants.vcf"

    add(haplotypeCaller)
  }
}
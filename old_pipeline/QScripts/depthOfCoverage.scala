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

import org.broadinstitute.sting.queue.extensions.snpeff._
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.apache.commons.io.FilenameUtils
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion


class depthOfCov extends QScript {
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R", required=true)
  var referenceFile: File = _ // _ is scala shorthand for null
  //var referenceFile: File = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta" // testing ref genome

  @Input(doc="One or more bam files.", shortName="I")
  var bamFiles: List[File] = Nil

  @Input(doc="Output core filename.", shortName="o", required=true)
  var out: File = _

  // The following arguments are all optional.
  @Input(doc="Calculate coverage statistics over this list of genes", shortName="geneList", required=false)
  var geneList: File = _
  
  @Argument(doc="Minimum mapping quality of reads to count towards depth.", shortName="mmq", required=false)
  var mmq: Int = -1
  
  @Argument(doc="for summary file outputs, report the % of bases coverd to >= this number.", shortName="ct", required=false)
  var ct: Seq[String] = Seq("15")

  @Input(doc="One or more genomic intervals over which to operate.", shortName="L", required=false)
  var intervalList: File = _

  // This trait allows us set the variables below in one place and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.memoryLimit = 4
  }

  def script() {
    // Create functions that we may run depending on options.
    val depthOfCoverage = new DepthOfCoverage with UnifiedGenotyperArguments

    // Create target list of intervals to be realigned
    // depthOfCoverage.num_threads = 12 -> only works with omitIntervalStatistics = True
    depthOfCoverage.input_file = bamFiles
    depthOfCoverage.R = referenceFile
    depthOfCoverage.calculateCoverageOverGenes = geneList
    depthOfCoverage.minMappingQuality = mmq
    depthOfCoverage.summaryCoverageThreshold = ct
    depthOfCoverage.intervals :+= intervalList
    depthOfCoverage.out = qscript.out + ".txt"

    add(depthOfCoverage)
  }
}
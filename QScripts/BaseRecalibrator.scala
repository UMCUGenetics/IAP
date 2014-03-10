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

class Realigner extends QScript {
    // Create an alias 'qscript' to be able to access variables in the Realigner.
    // 'qscript' is now the same as 'Realigner.this'
    qscript =>

    // Required arguments.  All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One bam file.", shortName="I",required=true)
    var bamFile: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    // Optional arguments
    @Input(doc ="Database of known polymorphic sites to skip over in the recalibration algorithm", shortName="knownSites", required=false)
    var knownFiles: List[File] = Nil
    
    // This trait allows us set the variables below in one place,
    // and then reuse this trait on each CommandLineGATK function below.
    trait BR_Arguments extends CommandLineGATK {
	this.reference_sequence = qscript.referenceFile
	this.memoryLimit = maxMem
    }

    def script() {
	val baseRecalibrator = new BaseRecalibrator with BR_Arguments
	val baseRecalibratorSecond = new BaseRecalibrator with BR_Arguments
	val analyzeCovariates = new AnalyzeCovariates with BR_Arguments
	val printReads = new PrintReads with BR_Arguments

	// Analyze patterns of covariation in the sequence dataset
	baseRecalibrator.input_file :+= bamFile
	if(knownFiles != Nil){
	    baseRecalibrator.knownSites = knownFiles
	}
	baseRecalibrator.out = swapExt(bamFile, ".bam", "_recal_data.table")
	
	baseRecalibrator.scatterCount = numScatters
	baseRecalibrator.nct = numCPUThreads
	
	// Do a second pass to analyze covariation remaining after recalibration
	baseRecalibratorSecond.input_file :+= bamFile
	if(knownFiles != Nil){
	    baseRecalibratorSecond.knownSites = knownFiles
	}
	baseRecalibratorSecond.BQSR = baseRecalibrator.out
	baseRecalibratorSecond.out = swapExt(bamFile, ".bam", "_post_recal_data.table")
	
	baseRecalibratorSecond.scatterCount = numScatters
	baseRecalibratorSecond.nct = numCPUThreads

	// Generate before and after plots
	analyzeCovariates.before = baseRecalibrator.out
	analyzeCovariates.after = baseRecalibratorSecond.out
	analyzeCovariates.plots = swapExt(baseRecalibrator.out, "recal_data.table", "baseRecalibration.pdf")

	// Apply the recalibration to your sequence data
	printReads.input_file :+= bamFile
	printReads.BQSR = baseRecalibrator.out
	printReads.out = swapExt(bamFile, ".bam", "_recalibrated.bam")

	printReads.scatterCount = numScatters
	printReads.nct = numCPUThreads

	add(baseRecalibrator,baseRecalibratorSecond,analyzeCovariates,printReads)
    }
}
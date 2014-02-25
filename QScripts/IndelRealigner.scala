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
import org.broadinstitute.sting.queue.function.ListWriterFunction

class Realigner extends QScript {
    // Create an alias 'qscript' to be able to access variables in the Realigner.
    // 'qscript' is now the same as 'Realigner.this'
    qscript =>

    // Required arguments.  All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _ 

    @Input(doc="One or more bam files.", shortName="I",required=true)
    var bamFiles: List[File] = Nil

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nt", required=true)
    var numDataThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _
    
    @Argument(doc="Realign mode: single or multi sample.", shortName="mode", required=true)
    var realignMode: String = _
    
    // Optional arguments
    @Input(doc ="Input VCF file with known indels.", shortName="known", required=false)
    var knownIndelFiles: List[File] = Nil
    
    // This trait allows us set the variables below in one place,
    // and then reuse this trait on each CommandLineGATK function below.
    trait TCIR_Arguments extends CommandLineGATK {
	this.reference_sequence = qscript.referenceFile
	this.memoryLimit = maxMem
    }

    def script() {
	if (realignMode == "multi") {
	    val targetCreator = new RealignerTargetCreator with TCIR_Arguments
	    val indelRealigner = new IndelRealigner with TCIR_Arguments

	    //Target creator
	    targetCreator.input_file = bamFiles
	    if(knownIndelFiles != Nil){
		targetCreator.known = knownIndelFiles
	    }
	    targetCreator.num_threads = numDataThreads
	    targetCreator.scatterCount = numScatters
	    targetCreator.out = new File("target.intervals.list")

	    //Indel realigner
	    indelRealigner.targetIntervals = targetCreator.out
	    indelRealigner.input_file = bamFiles
	    indelRealigner.scatterCount = numScatters
	    indelRealigner.nWayOut = "_realigned.bam"
	    add(targetCreator,indelRealigner)
	
	} else if (realignMode == "single") {
	    for (bamFile <- bamFiles) {
		val targetCreator = new RealignerTargetCreator with TCIR_Arguments
		val indelRealigner = new IndelRealigner with TCIR_Arguments

		//Target Creator
		targetCreator.input_file :+= bamFile
		if(knownIndelFiles != Nil){
		    targetCreator.known = knownIndelFiles
		}
		targetCreator.num_threads = numDataThreads
		targetCreator.scatterCount = numScatters
		targetCreator.out = swapExt(bamFile, "bam", "target.intervals.list")
		
		//Indel realigner
		indelRealigner.targetIntervals = targetCreator.out
		indelRealigner.input_file +:= bamFile
		indelRealigner.scatterCount = numScatters
		indelRealigner.out = swapExt(bamFile, "bam", "realigned.bam")

		add(targetCreator, indelRealigner)
	    }
	}
    }
}
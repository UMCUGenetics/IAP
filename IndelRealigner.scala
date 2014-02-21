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

package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.ListWriterFunction

class Realigner extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ 

  @Input(doc="One or more bam files.", shortName="I")
  var bamFiles: List[File] = Nil

  @Argument(doc="Scattercount.", shortName="SC")
  var scatter: Int = _
  
  @Argument(doc="Maxmem.", shortName="MEM")
  var maxmem: Double = _
  
  @Argument(doc="Number of threads.", shortName="TH")
  var threads: Int = _
  

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait TC_Arguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = qscript.maxmem
    this.num_threads = qscript.threads
  }

  trait IR_Arguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = qscript.maxmem
  }

  def script() {
    val targetCreator = new RealignerTargetCreator with TC_Arguments
    val indelRealigner = new IndelRealigner with IR_Arguments
    

    targetCreator.scatterCount = qscript.scatter
    targetCreator.input_file = qscript.bamFiles
    targetCreator.out = new File("target_intervals.list")
    
    indelRealigner.scatterCount = qscript.scatter
    indelRealigner.input_file = qscript.bamFiles
    indelRealigner.targetIntervals = targetCreator.out
    indelRealigner.nWayOut = "_realigned.bam"
    

    add(targetCreator, indelRealigner)


  }
}
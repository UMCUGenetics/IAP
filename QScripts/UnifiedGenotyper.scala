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

class VariantCaller extends QScript {
    // Create an alias 'qscript' to be able to access variables in the VariantCaller.
    // 'qscript' is now the same as 'VariantCaller.this'
    qscript =>

    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One or more bam files.", shortName="I")
    var bamFiles: List[File] = Nil

    @Input(doc="Output core filename.", shortName="O", required=true)
    var out: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _
    
    @Argument(doc="Genotype likelihoods calculation model to employ (SNP, INDEL or BOTH)", shortName="glm", required=true)
    var glm: String = _

    // The following arguments are all optional.
    @Argument(doc="Minimum phred-scaled confidence to call variants", shortName="stand_call_conf", required=false)
    var standCallConf: Int = 30 //default: best-practices value

    @Argument(doc="Minimum phred-scaled confidence to emit variants", shortName="stand_emit_conf", required=false)
    var standEmitConf: Int = 10 //default: best-practices value

    @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
    var dbsnpFile: File = _

    @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
    var targetFile: File = _

    def script() {
	val unifiedGenotyper = new UnifiedGenotyper

	// All required input
	unifiedGenotyper.input_file = bamFiles
	unifiedGenotyper.reference_sequence = referenceFile

	unifiedGenotyper.scatterCount = numScatters
	unifiedGenotyper.memoryLimit = maxMem
	unifiedGenotyper.num_cpu_threads_per_data_thread = numCPUThreads
	
	//SNP INDEL or BOTH
	if(glm == "SNP"){
	    unifiedGenotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.SNP
	    unifiedGenotyper.out = qscript.out + ".UG.raw_SNP.vcf"
	} else if(glm == "INDEL"){
	    unifiedGenotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.INDEL
	    unifiedGenotyper.out = qscript.out + ".UG.raw_INDEL.vcf"
	} else if(glm == "BOTH"){
	    unifiedGenotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.BOTH
	    unifiedGenotyper.out = qscript.out + ".UG.raw_variants.vcf"
	}
	// Optional input
	unifiedGenotyper.stand_emit_conf = standEmitConf
	unifiedGenotyper.stand_call_conf = standCallConf

	if(dbsnpFile != null) {
	    unifiedGenotyper.D = dbsnpFile
	}
	if(targetFile != null) {
	    unifiedGenotyper.L :+= targetFile
	}

	//add function to queue
	add(unifiedGenotyper)
    }
}
Plot Illumina Metrics
====================
Generates pdf and html reports from picard metrics.

The tool assumes the following directory structure and files
====================
	+ mapped_run_dir
		+ sample_1
			- sample_1_HSMetrics.txt
			- sample_1_MultipleMetrics.txt.insert_size_metrics
			- sample_1_MultipleMetrics.txt.quality_by_cycle_metrics
			- sample_1_MultipleMetrics.txt.quality_distribution_metrics
		+ sample_2
			- sample_2_HSMetrics.txt
			- sample_2_MultipleMetrics.txt.insert_size_metrics
			- sample_2_MultipleMetrics.txt.quality_by_cycle_metrics
			- sample_2_MultipleMetrics.txt.quality_distribution_metrics

Usage
====================
    cd mapped_run_dir
    perl src_dir/plotIlluminaMetrics.pl
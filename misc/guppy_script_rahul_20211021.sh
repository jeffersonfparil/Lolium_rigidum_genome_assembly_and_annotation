~/ont-guppy/bin/guppy_basecaller \
	-c dna_r9.4.1_450bps_sup.cfg \
	-i fast5/ \
	-s fastq \
	-x 'auto' \
	--recursive \
	--chunks_per_runner 412 \
	--gpu_runners_per_device 18 \
	--cpu_threads_per_caller 18


for h in 0.5 1.0; do
	for t in 0.01 0.05 0.1 0.5; do
		srun \
		--nodes 1 \
		--ntasks 1 \
		--cpus-per-task 16 \
		--partition normal \
		--job-name AL \
		--output 'AL_out.txt' \
		--open-mode append \
			./allele_age_simulator \
				--population_size=1000 \
				--theta=$t \
				--dominance=$h \
				--observed=10 \
				--replicates=10000000 \
				--seed=1000 & 
	done
done

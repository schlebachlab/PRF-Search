#!/bin/bash
#SBATCH -A highmem
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --job-name AspN-CAD
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out

apptainer run /depot/bsdrown/apps/images/fragpipe_23.1-complete.sif \
	--headless \
	--workflow fragpipe-second-pass.workflow \
	--manifest fragpipe-files-second-pass.fp-manifest \
	--workdir AspN_CAD_2pass_prophet_results \
	--threads 32 \
	--ram 254

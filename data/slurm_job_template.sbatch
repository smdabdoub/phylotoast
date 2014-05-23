#SBATCH --job-name={job_name}_{job_num}

# Walltime Limit: hh:nn:ss 
#SBATCH --time={walltime} 

# Node Specification:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={threads}

# Send email notification upon job fail
#SBATCH --mail-type=FAIL

module load qiime
cp {database_path} $TMPDIR
cp {job_num}.fna $TMPDIR
cd $TMPDIR

/usr/bin/time parallel_pick_otus_blast.py -i {job_num}.fna -r {database_fname} -O {threads} -s {similarity} -o bpo.{job_num}

cp -r bpo.{job_num} $SLURM_SUBMIT_DIR

#!/bin/bash

#Submit this script with: sbatch filename

#SBATCH --time=06:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=22   # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=128   # number of tasks per node
#SBATCH --job-name=RIOcalcs   # job name
#SBATCH --account=t26_coastal_ocean   # account name
#SBATCH --qos=standard  # qos name
#SBATCH --mail-user=sprice@lanl.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue   # do not requeue when preempted and on node failure
#SBATCH --signal=23@60  # send signal to job at [seconds] before end

# LOAD MODULEFILES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source /users/sprice/miniforge3/etc/profile.d/conda.sh
conda activate /users/sprice/miniforge3/envs/e3sm-unified

SCRIPTS=(
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops0.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops1.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops3.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops4.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops5.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops6.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops7.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops8.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops9.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops10.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops11.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops12.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops13.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops14.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops15.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops16.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops17.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops18.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops19.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops20.py
  /lustre/scratch5/sprice/calcRIODailyStats_NewLoops21.py
)

# sanity check
if [ "${#SCRIPTS[@]}" -ne 21 ]; then
  echo "WARNING: SCRIPTS array length != 21 (found ${#SCRIPTS[@]})"
fi

# create logs dir
mkdir -p logs

# launch each script as its own exclusive srun task and background it
for script in "${SCRIPTS[@]}"; do
  echo "Starting $script"
  srun --exclusive -n1 --cpus-per-task=1 python "$script" \
       > "logs/$(basename $script).out" 2> "logs/$(basename $script).err" &
done

# wait for all background tasks to complete
wait

echo "All tasks finished."

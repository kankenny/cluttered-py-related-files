#!/bin/sh

#SBATCH --partition=defq
#SBATCH --ntasks=38

module add parallel shared gcc/7.2.0  ipopt/3.13.4

MY_PROGRAM=neutrinos_cpp


echo "`date` Starting job on host `hostname`"
echo "#############################################################################################################"
echo



# This specifies the options used to run srun. The "-N1 -n1" options are
# used to allocates a single core to each task.
srun="srun --exclusive -N1 -n1"

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog ${MY_PROGRAM}.joblog.log --resume"

# Run a script, runtask.sh, using GNU parallel and srun. Parallel
# will run the runtask script for the numbers 1 through 38. To
# illustrate, the first job will run like this:
#
#   srun --exclusive -N1 -n1 ./submit.sh arg1:1 > runtask.1
#
$parallel "$srun ./${MY_PROGRAM} {1} > ${MY_PROGRAM}.log.{1} 2>&1" ::: {1..4}

# Note that if your program does not take any input, use the -n0 option to
# call the parallel command:
#
#   $parallel -n0 "$srun ./run_noinput_task.sh > output.{1}" ::: {1..38}

echo
echo "#############################################################################################################"
echo "`date` Ending job on host `hostname`"


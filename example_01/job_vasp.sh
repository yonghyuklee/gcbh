#!/bin/sh -f
#$ -l h_rt=24:00:00,h_data=4G
#$ -cwd
#$ -o LOG.$JOB_NAME.$JOB_ID
#$ -j y
#$ -pe shared 32
#$ -m bea

export Nprocs=32
echo RUNNING ON $HOSTNAME
cat $PE_HOSTFILE

source /u/local/Modules/default/init/modules.sh
source ~/.bashrc
module purge
module load IDRE intel/2020.4 intel/mpi hdf5/1.12.0_intel2019.2 anaconda3/2020.11
conda activate cupd

export ASE_VASP_COMMAND="mpirun -n $Nprocs /u/home/y/ylee/code/vasp/vasp.6.4.1/vasp.6.4.1/bin/vasp_std"
export VASP_PP_PATH=/u/home/y/ylee/code/vasp/vasp.6.4.1/POTCAR

python run_bh_ase.py


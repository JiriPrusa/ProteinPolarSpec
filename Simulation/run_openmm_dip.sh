#!/bin/bash
#run like batch_openmm_gpu [name] [hours]

DIR=$(pwd)
dir=${PWD##*/}
NPROC=2
name=$1
hour=$2

if [ -e ${name}.bsh ]
then
rm ${name}.bsh
echo "previous bsh removed" 
fi
cat << END >> ${name}.bsh
#PBS -l select=1:ncpus=${NPROC}:ngpus=1:mem=8gb:scratch_local=8gb:cluster=^doom
#PBS -q gpu_long
#PBS -j oe
#PBS -l walltime=${hour}:00:00

export OMP_NUM_THREADS=\$PBS_NUM_PPN
export OPENMM_CUDA_COMPILER=/storage/plzen1/home/prusaj/scripts/conda_nvcc
trap "clean_scratch" TERM EXIT

module add conda-modules-py37
module add python36-modules-gcc
conda activate openMM
cd \$SCRATCHDIR || exit 1
cp -r $DIR . || exit 2
cd $dir

DATADIR="\$PBS_O_WORKDIR"

python ./3_production_meta.py -p jesslys_let.pdb -s state_number_${name}.xml -n 10000000 --dt 1 --nrespa 2 -ch 'A' -x 'traj.dcd' -d 'dip.csv' --interval 60 --checkpoint state.chk

cp -r * $DIR || export CLEAN_SCRATCH=false
END
qsub ${name}.bsh

#!/bin/bash

#
#SBATCH -D /nica/mpd1/parfenov/TMP/
#SBATCH -J mpdEPNoTPC_step2
#SBATCH --mem-per-cpu=2G
#SBATCH -p nica
#SBATCH --time=20:30:00
#SBATCH -a 1-1
#
#SBATCH -o /nica/mpd1/parfenov/TMP/slurm_mpdEPNoTPC_step2_%A_%a.out
#SBATCH -e /nica/mpd1/parfenov/TMP/slurm_mpdEPNoTPC_step2_%A_%a.err
#

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/
export MPDROOT=/lhep/users/parfenov/Soft/mpdroot/install
source /lhep/users/parfenov/Soft/mpdroot/install/config/env.sh

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

export EP_STEP1_FILE=/nica/mpd1/parfenov/mpdConverter/OUT/ep_step1_mpdtree_urqmd_xew_2.5agev_prod40/3944082/ep.root
export FILELIST=/lhep/users/parfenov/Soft/mpdConvert/macros/mpdtree_urqmd_xew_2.5agev_prod40.list

export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export LABEL=${SHORTNAME11}

export INFILE=`sed "${TASK_ID}q;d" ${FILELIST}`

if [[ -f "$INFILE" ]]; then
export DATE=${JOB_ID} # or `date '+%Y%m%d_%H%M%S'`

export MAIN_DIR=/lhep/users/parfenov/Soft/mpdConvert
export MAIN_OUT=/nica/mpd1/parfenov/mpdConverter
export OUT_DIR=${MAIN_OUT}/OUT/ep_step2_${LABEL}/${DATE}

export OUT_FILE_DIR=${OUT_DIR}/files
export OUT_LOG_DIR=${OUT_DIR}/log
export OUT_FILE=${OUT_FILE_DIR}/ep_step2_${LABEL}_${JOB_ID}_${TASK_ID}.root
export OUT_LOG=${OUT_LOG_DIR}/ep_step2_${LABEL}_${JOB_ID}_${TASK_ID}.log

export TMP=/nica/mpd1/parfenov/TMP
#${MAIN_DIR}/TMP
export TMP_DIR=${TMP}/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $OUT_FILE_DIR
mkdir -p $OUT_LOG_DIR
mkdir -p $TMP_DIR
touch $OUT_LOG

cp -v ${MAIN_DIR}/runEP_step2.c ${TMP_DIR}/runEP_step2.c &>> $OUT_LOG

echo "Input file : ${INFILE}" &>> $OUT_LOG
echo "Output file: ${OUT_FILE}" &>> $OUT_LOG

cd ${TMP_DIR}
time root -l -b -q runEP_step2.c'("'${INFILE}'","'${OUT_FILE}'","'${EP_STEP1_FILE}'")' &>> $OUT_LOG

rm -rfv ${TMP_DIR} &>> $OUT_LOG

echo "Job is finished!" &>> $OUT_LOG

fi

#!/bin/bash

#
#SBATCH -D /lustre/home/user/p/parfenov/TMP
#SBATCH -J mpdQA
#SBATCH --mem-per-cpu=2G
#SBATCH -p mephi
#SBATCH --time=02:30:00
#SBATCH -a 1-1
#
#SBATCH -o /lustre/home/user/p/parfenov/TMP/slurm_mpdQA_%A_%a.out
#SBATCH -e /lustre/home/user/p/parfenov/TMP/slurm_mpdQA_%A_%a.err
#

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v23.09.23-1
export MPDROOT=/lustre/home/user/p/parfenov/Soft/mpdroot/install
source /lustre/home/user/p/parfenov/Soft/mpdroot/install/config/env.sh

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

export ecm=3.5 #2.5, 3.0, 3.5 GeV
export nucl_mass=209 #209 for Bi+Bi

export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/mpdtree_urqmd_bibi_${ecm}gev_mpdfxt.list

export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export LABEL=${SHORTNAME11}

export INFILE=`sed "${TASK_ID}q;d" ${FILELIST}`

if [[ -f "$INFILE" ]]; then
export DATE=${JOB_ID} # or `date '+%Y%m%d_%H%M%S'`

export MAIN_DIR=/lustre/home/user/p/parfenov/Soft/mpdConvert
export OUT_DIR=${MAIN_DIR}/OUT/qa_${LABEL}/${DATE}
export OUT_FILE_DIR=${OUT_DIR}/files
export OUT_LOG_DIR=${OUT_DIR}/log
export OUT_FILE=${OUT_FILE_DIR}/qa_${LABEL}_${JOB_ID}_${TASK_ID}.root
export OUT_LOG=${OUT_LOG_DIR}/qa_${LABEL}_${JOB_ID}_${TASK_ID}.log

export TMP=${MAIN_DIR}/TMP
export TMP_DIR=${TMP}/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $OUT_FILE_DIR
mkdir -p $OUT_LOG_DIR
mkdir -p $TMP_DIR
touch $OUT_LOG

cp -v ${MAIN_DIR}/runQaMpd.C ${TMP_DIR}/runQaMpd.C &>> $OUT_LOG

echo "Input file : ${INFILE}" &>> $OUT_LOG
echo "Output file: ${OUT_FILE}" &>> $OUT_LOG

cd ${TMP_DIR}
time root -l -b -q runQaMpd.C'("'${INFILE}'","'${OUT_FILE}'","'${ecm}'", "'${nucl_mass}'")' &>> $OUT_LOG

rm -rfv ${TMP_DIR} &>> $OUT_LOG

echo "Job is finished!" &>> $OUT_LOG

fi

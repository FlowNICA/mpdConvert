#!/bin/bash

#
#SBATCH -D /lustre/stor2/mephi/parfenov/TMP
#SBATCH -J mpdQASt2
#SBATCH --mem-per-cpu=2G
#SBATCH -p mephi
#SBATCH --time=02:30:00
#SBATCH -a 1-1
#
#SBATCH -o /lustre/stor2/mephi/parfenov/TMP/slurm_mpdQA_%A_%a.out
#SBATCH -e /lustre/stor2/mephi/parfenov/TMP/slurm_mpdQA_%A_%a.err
#

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.12.24-1
export MPDROOT=/lustre/home/user/p/parfenov/Soft/mpdroot/install
source /lustre/home/user/p/parfenov/Soft/mpdroot/install/config/env.sh

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

export ecm=2.87 #2.5, 3.0, 3.5 GeV
export nucl1_mass=184 #209 for Bi, 184 for W, 154 for Xe
export nucl2_mass=154 #209 for Bi, 184 for W, 154 for Xe
export system=xew
export programm=prod36

export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/lists/mpdtree_urqmd_${system}_${ecm}gev_${programm}.list
#export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/mpdtree_urqmd_${system}_${ecm}gev_${programm}.list

#export PID_FILE=/lustre/stor2/mephi/parfenov/mpdtree/OUT/qa_step1_mpdtree_urqmd_xexe_2.87gev_prod36/8818308/pid.root
export PID_FILE=/lustre/stor2/mephi/parfenov/mpdtree/OUT/qa_step1_mpdtree_urqmd_xew_2.87gev_prod36/8818608/pid.root

export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export LABEL=${SHORTNAME11}

export INFILE=`sed "${TASK_ID}q;d" ${FILELIST}`

if [[ -f "$INFILE" ]]; then
export DATE=${JOB_ID} # or `date '+%Y%m%d_%H%M%S'`

export MAIN_DIR=/lustre/home/user/p/parfenov/Soft/mpdConvert
export OUT_DIR=/lustre/stor2/mephi/parfenov/mpdtree/OUT/qa_step2_${LABEL}/${DATE}
#${MAIN_DIR}/OUT/qa_${LABEL}/${DATE}
export OUT_FILE_DIR=${OUT_DIR}/files
export OUT_LOG_DIR=${OUT_DIR}/log
export OUT_FILE=${OUT_FILE_DIR}/qa_${LABEL}_${JOB_ID}_${TASK_ID}.root
export OUT_LOG=${OUT_LOG_DIR}/qa_${LABEL}_${JOB_ID}_${TASK_ID}.log

export TMP=/lustre/stor2/mephi/parfenov/TMP
#${MAIN_DIR}/TMP
export TMP_DIR=${TMP}/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $OUT_FILE_DIR
mkdir -p $OUT_LOG_DIR
mkdir -p $TMP_DIR
touch $OUT_LOG

cp -v ${MAIN_DIR}/runQaMpd_step2.C ${TMP_DIR}/runQaMpd_step2.C &>> $OUT_LOG
cp -v ${MAIN_DIR}/FunctionsQa.C ${TMP_DIR}/FunctionsQa.C &>> $OUT_LOG

echo "Input file : ${INFILE}" &>> $OUT_LOG
echo "Output file: ${OUT_FILE}" &>> $OUT_LOG
echo "PID file:    ${PID_FILE}" &>> $OUT_LOG

cd ${TMP_DIR}
time root -l -b -q runQaMpd_step2.C'("'${INFILE}'","'${OUT_FILE}'","'${ecm}'", "'${nucl1_mass}'", "'${nucl2_mass}'", "'${PID_FILE}'")' &>> $OUT_LOG

rm -rfv ${TMP_DIR} &>> $OUT_LOG

echo "Job is finished!" &>> $OUT_LOG

fi

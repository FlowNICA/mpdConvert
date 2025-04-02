#!/bin/bash

#
#SBATCH -D /lustre/stor2/mephi/parfenov/TMP/
#SBATCH -J ReConvertMpd
#SBATCH --mem-per-cpu=4G
#SBATCH -p cascade
#SBATCH --time=36:30:00
#SBATCH -a 1-1
#
#SBATCH -o /lustre/stor2/mephi/parfenov/TMP/slurm_mpdreconv_%A_%a.out
#SBATCH -e /lustre/stor2/mephi/parfenov/TMP/slurm_mpdreconv_%A_%a.out
#

ls /cvmfs/nica.jinr.ru/
sleep 5

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.12.24-1
export MPDROOT=/lustre/home/user/p/parfenov/Soft/mpdroot/install
source /lustre/home/user/p/parfenov/Soft/mpdroot/install/config/env.sh

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}


#export ecm=3.5
#export prefix=mpdfxt
export ecm=9.2
export prefix=req25

#export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/urqmd_xexe_${ecm}gev_${prefix}.list
#export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/urqmd_xew_${ecm}gev_${prefix}.list
export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/urqmd_bibi_${ecm}gev_${prefix}.list
#export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/phsd_bibi_${ecm}gev.list
#export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/prod36_urqmd_xew_2.87gev_mpdfxt.list
#export FILELIST=/lustre/home/user/p/parfenov/Soft/mpdConvert/macros/prod36_urqmd_xexe_2.87gev_mpdfxt.list

export GEOFILE=${MPDROOT}/geometry/zdc_oldnames_7sect_v1_no_overlaps_w_pipe_magnet.root

export ORIGINAL_OUT_DIR=$1

export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export LABEL=${SHORTNAME11}

export INFILE=`sed "${TASK_ID}q;d" ${FILELIST}`

if [[ -f "$INFILE" ]]; then
export ORIGINAL_JOB_ID=`echo ${ORIGINAL_OUT_DIR} | cut -d"/" -f8`
export orig2=`echo $orig1 | cut -d"/" -f8`
export ORIGINAL_OUT_FILE=${ORIGINAL_OUT_DIR}/files/mpdtree_${ORIGINAL_JOB_ID}_${TASK_ID}.root
export ORIGINAL_OUT_LOG=${ORIGINAL_OUT_DIR}/log/mpdtree_${ORIGINAL_JOB_ID}_${TASK_ID}.log

if [[ ! -f "$ORIGINAL_OUT_FILE" ]]; then

export MAIN_DIR=/lustre/home/user/p/parfenov/Soft/mpdConvert
export MAIN_OUT=/lustre/stor2/mephi/parfenov/mpdtree
export OUT_DIR=${MAIN_OUT}/OUT/${LABEL}/${DATE}

export OUT_FILE_DIR=${OUT_DIR}/files
export OUT_LOG_DIR=${OUT_DIR}/log
export OUT_FILE=${OUT_FILE_DIR}/mpdtree_${LABEL}_${ORIGINAL_JOB_ID}_${TASK_ID}.root
export OUT_LOG=${OUT_LOG_DIR}/mpdtree_${LABEL}_${ORIGINAL_JOB_ID}_${TASK_ID}.log

export TMP=/lustre/stor2/mephi/parfenov/TMP
export TMP_DIR=${TMP}/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $OUT_FILE_DIR
mkdir -p $OUT_LOG_DIR
mkdir -p $TMP_DIR
#touch $OUT_LOG

cp -v ${MAIN_DIR}/convertMPD.C ${TMP_DIR}/convertMPD.C &>> $OUT_LOG

source /scratch2/parfenov/Soft/mpdroot/install/config/env.sh

echo "" &>> $OUT_LOG
echo "--------------------Restarting job!-------------------" &>> $OUT_LOG
echo "" &>> $OUT_LOG
echo "Input file : ${INFILE}" &>> $OUT_LOG
echo "Original Job_Id: ${ORIGINAL_JOB_ID}" &>> $OUT_LOG
echo "Output file: ${OUT_FILE}" &>> $OUT_LOG

#root -l -b -q ${TMP_DIR}/convertMPD.C+'("'${TMP_DIR}/input.list'","'$OUT_FILE'")' &>> $OUT_LOG
#root -l -b -q ${TMP_DIR}/convertMPD.C+'("'${TMP_DIR}/mpddst.root'","'$OUT_FILE'")' &>> $OUT_LOG
root -l -b -q ${TMP_DIR}/convertMPD.C+'("'${INFILE}'","'$OUT_FILE'")' &>> $OUT_LOG

rm -rfv ${TMP_DIR} &>> $OUT_LOG

echo "Job is finished!" &>> $OUT_LOG

fi
fi

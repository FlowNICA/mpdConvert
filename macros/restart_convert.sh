#!/bin/bash

#
#$ -wd /scratch2/$USER/TMP
#$ -cwd
#$ -l h=!(ncx152.jinr.ru|ncx205.jinr.ru|ncx123.jinr.ru|ncx111.jinr.ru|ncx113.jinr.ru|ncx149.jinr.ru|ncx223.jinr.ru)
#$ -N rerun_mpddst2tree
#$ -q all.q
#$ -l h_rt=0:10:00
#$ -l s_rt=0:10:00
#$ -t 1-20000
#
#$ -o /scratch2/$USER/TMP
#$ -e /scratch2/$USER/TMP
#

source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev

export JOB_ID=${JOB_ID}
export TASK_ID=${SGE_TASK_ID}

export ecm=3.5

export FILELIST=/scratch2/parfenov/Soft/mpdConvert/macros/lists/urqmd_bibi_${ecm}gev_mpdfxt.list
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

export MAIN_DIR=/scratch2/parfenov/Soft/mpdConvert
export OUT_DIR=${ORIGINAL_OUT_DIR}
export OUT_FILE_DIR=${OUT_DIR}/files
export OUT_LOG_DIR=${OUT_DIR}/log
export OUT_FILE=${OUT_FILE_DIR}/mpdtree_${LABEL}_${ORIGINAL_JOB_ID}_${TASK_ID}.root
export OUT_LOG=${OUT_LOG_DIR}/mpdtree_${LABEL}_${ORIGINAL_JOB_ID}_${TASK_ID}.log

export TMP=${MAIN_DIR}/TMP
export TMP_DIR=${TMP}/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $OUT_FILE_DIR
mkdir -p $OUT_LOG_DIR
mkdir -p $TMP_DIR
#touch $OUT_LOG

#xrdcp $INFILE ${TMP_DIR}/mpddst.root &>> $OUT_LOG
rsync -vuzhP ${MAIN_DIR}/convertMPD.C ${TMP_DIR}/convertMPD.C &>> $OUT_LOG

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

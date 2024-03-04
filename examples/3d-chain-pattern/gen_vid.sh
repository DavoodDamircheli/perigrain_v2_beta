#!/bin/bash
OUTPUT_DIR=${1}
# scale while making video from pic
scale=1
IMG_PREF='vp_'

RNDN=$(date +%F-%H-%M-%S)-$(cat /dev/urandom | tr -cd 'a-z0-9' | head -c 5)

FILE=${OUTPUT_DIR}/${RNDN}.mp4

# standard
# framerate=20
# r=25
# crf=30

framerate=40
r=25
crf=30

ffmpeg -framerate $framerate -i ${OUTPUT_DIR}/${IMG_PREF}%5d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r $r -crf $crf ${FILE}

mplayer $FILE

echo $FILE

#!/bin/bash
OUTPUT_DIR=${1}
n=${2}  # This is the variable n to control how many images to process

# scale while making video from pic
scale=0.5
IMG_PREF='img_'

RNDN=$(date +%F-%H-%M-%S)-$(cat /dev/urandom | tr -cd 'a-z0-9' | head -c 5)

FILE=${OUTPUT_DIR}/${RNDN}.mp4

# standard
# framerate=20
# r=25
# crf=30

framerate=40
r=25
crf=30

# Pass the variable n to limit the number of frames
ffmpeg -framerate $framerate -i ${OUTPUT_DIR}/${IMG_PREF}%5d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -frames:v $n -r $r -crf $crf ${FILE}

mplayer $FILE

echo $FILE


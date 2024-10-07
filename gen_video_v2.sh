#!/bin/bash                                                                                 
dir=${1}                                                                                    
output_dir=${2}                                                                             
n1=${3}  # Starting image number (e.g., 1)                                                  
n2=${4}  # Ending image number (e.g., 100)                                                  
                                                                                            
# scale while making video from pic                                                         
scale=0.5                                                                                   
IMG_PREF='img_'                                                                             

RNDN=$(date +%F-%H-%M-%S)-$(cat /dev/urandom | tr -cd 'a-z0-9' | head -c 5)                 

FILE=${output_dir}/${RNDN}.mp4                                                              

# Calculate the number of frames (n2 - n1 + 1)                                              
frames=$(($n2 - $n1 + 1))                                                                   

# standard                                                                                  
framerate=40                                                                                
r=25                                                                                        
crf=30                                                                                      

# Create video from images starting at n1 and ending at n2                                  
ffmpeg -framerate $framerate -start_number $n1 -i ${dir}/${IMG_PREF}%5d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -frames:v $frames -r $r -crf $crf ${FILE}                                    

mplayer $FILE                                                                               

echo $FILE


#!/bin/bash
set eu
cd ~
 

ymdhms=$(date "+%Y-%m-%d_%H-%M-%S")
ymdhmsforcalc=$(date "+%s")
unixtime=$(date "+%s")
 
target_dir=$HOME/Documents/qfa/QFA${unixtime}_${ymdhms}
 

mkdir -p $target_dir
cd $target_dir

counter=1
while :
do
 counter_p=$(printf "%03d" $counter)
 echo $counter_p time epson-lower
 nowtime=$(date "+%Y-%m-%d_%H-%M-%S")
 
 # mirrorはnoにしているのは、Pythonでの解析時にflipするから
 scanimage --device "epson2:libusb:001:005" --resolution 600 --mode Gray --format=tiff --mirror=no > QFA${unixtime}_${nowtime}.tiff &

 counter=$(($counter + 1))
 sleep 10m
done


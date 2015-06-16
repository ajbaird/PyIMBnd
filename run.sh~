#! /bin/sh

# a small shell script to clean up the files in the folder and run the 
# immersed boundary code


echo 'removing folders and .png files'

rm *.png
rm -r velu velv vorticity forcesx forcesy pressure

echo 'running code'
python ibcode.py

#cp -r velu,velv,vorticity,forcesx,forcesy,pressure ~/Desktop/working_directory 



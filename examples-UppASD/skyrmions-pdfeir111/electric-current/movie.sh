#!/bin/bash
#Programer: Ivan Miranda

fac=1

rm -f current.mp4 2>/dev/null

lin_ncell=$(grep -in "ncell" inpsd.dat | awk -F: '{print $1}');
dim=$(sed "${lin_ncell}"'!d' inpsd.dat | awk '{print $2}');
num_at=$(echo "${dim}*${dim}" | bc);
num_lin=$(wc -l moment*.out | awk '{print $1}')
num_lin=$((${num_lin}-7))
num_steps=$(echo "(${num_lin})/(${fac}*${num_at})" | bc)

echo "Numero de atomos: ${num_at}"
echo "Numero de steps: ${num_steps}"
num_dig=$(echo "${#num_steps}");

for i in $(seq 0 ${num_steps}); do 

    echo -ne "Executando imagens... ${i}/${num_steps} \033[0K\r"

    if [ $(echo "${i} == 0" | bc) -eq 1 ] ; then
        k=8 ; j=$((${num_at}+7))
    else
        k=$(((${i}*${fac}-1)*${num_at}+8)); 
        j=$((${i}*${num_at}*${fac}+7)); 
    fi 

    sed "${k},${j}"'!d' moment*.out > restart.ni${i} ; 
    text=$(grep "restart*" ASD_snapmap.py | cut -d'"' -f 2 | tail -1); 
    sed -i 's/'"${text}"'/'"restart.ni${i}"'/g' ASD_snapmap.py ; ./ASD_snapmap.py ; 
    foo=$(printf "%0${num_dig}d" ${i});
    mv snapmap*${i}*.png ${foo}.png
    rm -f restart.ni${i} 2>/dev/null

done 

sed -i 's/'"restart.ni${num_steps}"'/restart.PdFeIr11.out/g' ASD_snapmap.py
cp 00.png 01.png ; cp 10.png 11.png
ffmpeg -r 15 -i %0${num_dig}d.png -c:v libx264 -r 30 -pix_fmt yuv420p current.mp4

echo
echo "Done!"


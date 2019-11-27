#!/bin/sh

./create_network -n 3000
./d_opf -s > iterates_cold.txt

END=100
for i in `seq 1 $END`
do
    scale=$(echo "$i*0.01" | bc)
    
    cp network network_save
    ./perturb_network -s $scale -f network_save

    filename=$(printf "iterates/iterates_warm_%03d.txt" $i)
    
    ./d_opf -w -f network_save > $filename
done

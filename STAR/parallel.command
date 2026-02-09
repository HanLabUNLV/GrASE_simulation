#parallel -j 12 bash STAR/02_pass1.sh group1 {} :::: group.list
#parallel -j 12 bash STAR/02_pass1.sh group2 {} :::: group.list
parallel -j 12 bash STAR/08_pass2.sh group1 {} :::: group.list
parallel -j 12 bash STAR/08_pass2.sh group2 {} :::: group.list

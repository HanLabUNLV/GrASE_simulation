
#cat *.final.SJ.tab \
#  | awk '{print $1"\t"$2"\t"$3}' \
#  | sort -k1,1 -k2,2n -k3,3n \
#  | uniq > all_celltypes.final.merged.SJ.tab1


cat *.final.SJ.tab \
  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
  | sort -k1,1 -k2,2n -k3,3n \
  | uniq > all_celltypes.final.merged.SJ.tab


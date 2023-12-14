## ref point 
echo 102.9 37.2 | proj +proj=utm +zone=48
313629.46   4119124.88

# create the station catalog with name and position
grep -v "#" Supp_Table_S1.txt | awk '{print $2, $3, $1}'|proj +proj=utm +zone=48 | awk '{printf("%s %f %f\n",$3,($1-313629.46)/ 1e3,($2-4119124.88)/1e3)}' > stations_liang_km.dat

# 
grep -v "#" table_liang_sblock_ll.dat | awk '{fname="./sblock/"$1;print 1,$4,$5,$6,$7,$8 > fname;close(fname);}'

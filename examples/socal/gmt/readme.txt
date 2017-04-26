# create the coastline files

pscoast -R-116.5/-115.5/33/33.5 -W -Jm -m -Df >ca_coast_ll.dat
cat ca_coast_ll.dat | awk '{if ($1 == ">") print ">"; else printf("%f %f\n",$1,$2)}' |proj +proj=utm +zone=11| awk '{if  ($1 == "*") print ">"; else printf("%f %f\n",($1-630639.72)/ 1e3,($2-3663244.76)/1e3, 0)}'> ca_coasts_km.xyz

# reduction gmt file
awk '{if ($1 == ">") print ">"; else {if ($1>-300 && $1<300 && $2>-300 && $2<300) {print $1, $2, $3}}}' ca_faults_km.xyz > ca_faults_300km.xyz

# conversion moho.xyz to csv
awk 'BEGIN{print "\"East\",\"North\",\"Down\",\"Depth\""}{print $1","$2","$3","$3}' cal_moho01_q02_q08_ir02_id01_km.xyz > cal_moho01_q02_q08_ir02_id01_km.csv

#conversion lat lon to km
grep -v "#" srtm_socal.xyz | proj +proj=utm +zone=11 | awk '{printf("%f %f  %f\n",($1-  630639.72)/1e3, ($2-3663244.76)/1e3 , $3)}' > srtm_socal_km.xyz

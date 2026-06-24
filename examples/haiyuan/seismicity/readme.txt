# convert csv from usgs to xydm
awk -F "," '{print $3,$2,$4,$5}' usgs_2000-2014_101-105_34-39.csv > usgs_2000-2014_101-105_34-39.xydm

## ref point 
echo 102.9 37.2 | proj +proj=utm +zone=48
313629.46   4119124.88

#convert lat lon to km
grep -v "#" usgs_2000-2014_101-105_34-39.xydm | proj +proj=utm +zone=48 | awk '{printf("%f %f %f\n"($1-313629.46)/1e3,($2-4119124.88)/1e3,$3,$4)}' > usgs_2000-2014_101-105_34-39_km.xydm

# convert cgps to list of stations
grep -v "#" pbo.final_frame.psvelo | awk '{print $1,$2,$8}' | proj +proj=utm +zone=11 | awk '{printf("%s %f %f\n",$3,($1-630639.72)/1e3,($2-3663244.76)/1e3)}' > cgps_stations_km.dat

# create a list of GPS positions for Paraview
ifile=cgps_stations_km.dat;ofile=cgps_stations_km.csv;echo "\"North(km)\",\"East(km)\",\"Depth\",\"NAME\"" > $ofile; awk '{print $3","$2",0,"$1}' $ifile >> $ofile

# create the data file for each station
grep -v "#" pbo.final_frame.psvelo | awk '{fname="./cgps2_interseismic/"$8;print $3,$4,$5,$6 > fname;close(fname);}'

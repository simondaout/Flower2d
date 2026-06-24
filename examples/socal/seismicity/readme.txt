# Hauksson+12
cat sc_02_11_geo.xydma | proj +proj=utm +zone=11 +ellps=WGS84 | awk '{printf("%f %f %f %f %d\n",($1-630639.72)/1e3,($2-3663244.76)/1e3,$3,$4,$5)}' > sc_02_11_km.xydm

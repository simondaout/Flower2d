# sampling
grdsample T104_spectrum-square_clean2_noflata_LOSVelocity_myr.grd -I200+/200+ -GT104_spectrum-square_clean_noflata_LOSVelocity_myr_s200.grd

# convert grd to xyz
grd2xyz T104_spectrum-square_clean_noflata_LOSVelocity_myr_s200.grd -S > T104_spectrum-square_clean_noflata_LOSVelocity_myr_s200.xylos

#lon lat los to km km los for insar file
grep -v "#" T104_spectrum-square_clean2_noflata_LOSVelocity_myr_s200.xylos | proj +proj=utm +zone=48 | awk '{printf("%f %f %f\n",($1-313629.46)/1e3, ($2-4119124.88)/1e3 , $3*-1000)}' > T104_spectrum-square_clean2_noflata_LOSVelocity_mmyr_s200_km.xy-los




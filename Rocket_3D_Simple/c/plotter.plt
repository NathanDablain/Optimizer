set term png
set output "Position_NED.png"
set xlabel "Time (s)"
set multiplot layout 3,1 rows
set ylabel "Position - North (m)"
plot 'Rocket_3D_Simple_log.txt' using 1:2 with lines linewidth 3 title ""
set ylabel "Position - East (m)"
plot 'Rocket_3D_Simple_log.txt' using 1:3 with lines linewidth 3 title ""
set ylabel "Position - Up (m)"
plot 'Rocket_3D_Simple_log.txt' using 1:(-$4) with lines linewidth 3 title ""
unset multiplot

set output "Speed_c_L.png"
set xlabel "Time (s)"
set multiplot layout 3,1 rows
set ylabel "Speed (m/s)"
plot 'Rocket_3D_Simple_log.txt' using 1:5 with lines linewidth 3 title ""
set ylabel "c_Ly"
plot 'Rocket_3D_Simple_log.txt' using 1:12 with lines linewidth 3 title ""
set ylabel "c_Lz"
plot 'Rocket_3D_Simple_log.txt' using 1:13 with lines linewidth 3 title ""
unset multiplot

unset term
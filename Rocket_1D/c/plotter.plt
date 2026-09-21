set term png
set output "h_v_T.png"
set xlabel "Time (s)"
set multiplot layout 3,1 rows
set ylabel "Height (m)"
plot 'Rocket1D_log.txt' using 1:2 with lines linewidth 3 title ""
set ylabel "Velocity (m/s)"
plot 'Rocket1D_log.txt' using 1:3 with lines linewidth 3 title ""
set ylabel "Thrust (N)"
plot 'Rocket1D_log.txt' using 1:5 with lines linewidth 3 title ""
unset multiplot
unset term
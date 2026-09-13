set term png
set output "p_v_a.png"
set xlabel "Time (s)"
set multiplot layout 3,1 rows
set ylabel "Position (m)"
plot 'Block1D_log.txt' using 1:2 with lines linewidth 3 title ""
set ylabel "Velocity (m/s)"
plot 'Block1D_log.txt' using 1:3 with lines linewidth 3 title ""
set ylabel "Acceleration (m/s/s)"
plot 'Block1D_log.txt' using 1:4 with lines linewidth 3 title ""
unset multiplot
unset term
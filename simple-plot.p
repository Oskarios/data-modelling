set xlabel "Iterations"
set ylabel "Model Parameters C and m"
set title "Evolution of model parameters for small dataset with greed parameter 0.1 and tolerance e-7"

# set yrange [-10:250]
# set xrange [0:300]
# set xscale

set yrange [2:10]
set xrange [-10:1400]

unset logscale
data = 'data_out/simple_set_out.dat'

# p data using 1:5, data using 1:6

p data using 1:2 w l title "C (initial: 10.0)", data using 1:3 w l title "m (initial: 10.0)", 3.22199,2.82199,9.58307,7.58307 
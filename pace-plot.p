set xlabel "Iterations"
set ylabel "Model Parameters C and m"
set title "Evolution of model parameters with greed parameter 0.01 and tolerance e-6"

set xrange [-5:7550]
set yrange [-0.7:8.2]

data = "data_out/pace_set_out.dat"

p data using 1:2 w l title "m (initial: 10.0)", data using 1:3 w l title "b (initial: 10.0)",0.860936,0.0409809
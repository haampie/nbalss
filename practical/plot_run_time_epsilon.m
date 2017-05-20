function plot_run_time_epsilon()
  epss = [1e-08, 1e-07, 5e-07, 1e-06, 2.5e-06, 5e-06];
  y = [64.83, 52.40, 41.56, 40.06, 71.15, 89.83];

  semilogx(epss, y, '-.O')
  grid on
  xlabel('\epsilon')
  ylabel('Run time (s)')

end
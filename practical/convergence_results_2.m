function convergence_results_2
  p = 2;
  x = [32.0, 64.0, 96.0, 128.0, 160.0];
  inf_norm = [1.02851, 1.07422, 1.08176, 1.08405, 1.08572];
  een_achtste = [0.327567, 0.381462, 0.401002, 0.410911, 0.416887];
  drie_achtste = [0.500404, 0.489258, 0.481327, 0.476320, 0.472971];



  inf_interp = richardson(x, inf_norm, p)
  err_inf_norm = abs(inf_norm - inf_interp);
  rel_inf = abs(inf_norm - inf_interp) ./ inf_norm

  second_order = 1 ./ x .^ p;
  subplot(1, 3, 1)
  loglog(x, err_inf_norm, '-*'); hold on;
  loglog(x, log_schaling(err_inf_norm, second_order)); hold off;
  legend('Approximate error', sprintf('~ 1 / N^%d', p))
  xlabel('N')
  grid on;
  title('Global (inf norm)')

  een_achtste_interp = richardson(x, een_achtste, p)
  err_een_achtste = abs(een_achtste - richardson(x, een_achtste, p));
  rel_een_achtste = abs(een_achtste - een_achtste_interp) ./ een_achtste
  second_order = 1 ./ x .^ p;
  subplot(1, 3, 2)
  loglog(x, err_een_achtste, '-*'); hold on;
  loglog(x, log_schaling(err_een_achtste, second_order)); hold off;
  legend('Approximate error', sprintf('~ 1 / N^%d', p))
  xlabel('N')
  grid on;
  title('Local (1/8)')

  drie_achtste_interp = richardson(x, drie_achtste, p)
  err_drie_achtste = abs(drie_achtste - drie_achtste_interp);
  rel_drie_achtste = abs(drie_achtste - drie_achtste_interp) ./ drie_achtste
  second_order = 1 ./ x .^ p;
  subplot(1, 3, 3)
  loglog(x, err_drie_achtste, '-*'); hold on;
  loglog(x, log_schaling(err_drie_achtste, second_order)); hold off;
  legend('Approximate error', sprintf('~ 1 / N^%d', p))
  xlabel('N')
  grid on;
  title('Local (3/8)')

end

function v = richardson(x, y, p)
  r = x(end) / x(end - 1);
  v = (r^p * y(end) - y(end - 1)) / (r^p - 1);
end

function s = log_schaling(y1, y2)
  s = (y1(end - 2) / y2(end - 2)) * y2;
end

function convergence_results
  % Global
  x = [32.0, 64.0, 96.0, 128.0, 160.0];
  inf_norm = [1.02851, 1.07422, 1.08176, 1.08405, 1.08572];

  % Local
  een_achtste = [0.327567, 0.381462, 0.401002, 0.410911, 0.416887];
  drie_achtste = [0.500404, 0.489258, 0.481327, 0.476320, 0.472971];

  figure;
  diff1 = difference1(inf_norm);
  loglog(x(1 : end - 1), diff1, 'b-*')
  order_inf_norm = the_order(x, diff1)
  title('Global (inf norm)')

  figure;
  diff2 = difference1(een_achtste);
  loglog(x(1 : end - 1), diff2, 'b-*')
  order_een_achtste = the_order(x, diff2)
  title('Local (1/8)')

  figure;
  diff3 = difference1(drie_achtste);
  loglog(x(1 : end - 1), diff3, 'b-*')
  order_drie_achtste = the_order(x, diff3)
  title('Local (3/8)')

end

function d = difference1(vec)
  d = abs(vec(2 : end) - vec(1 : end - 1));
end

function d = difference2(vec)
  d = abs(vec(1 : end - 1) - vec(end));
end

function r = the_order(xs, ys)
   r = abs((log(ys(2 : end)) - log(ys(1 : end - 1))) ./ (log(xs(2 : end - 1)) - log(xs(1 : end - 2))));
end
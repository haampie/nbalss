function cont(mu, theta, mu_end, steps, tol)
  mus = linspace(mu, mu_end, steps);

  legend_entries = {};
  
  figure;
  hold on;
  for mu = mus
    [theta, history] = zeroOfBuck(length(theta), mu, tol, theta);
    fprintf('Steps = %d\n', length(history));
    legend_entries{end + 1} = ['Mu = ' num2str(mu)];
    plot(theta);
  end
  legend(legend_entries)
end
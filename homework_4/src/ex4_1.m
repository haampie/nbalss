function ex4_1
  n = 800;
  thetas = init(10, n, 1e-6);
  at_mu_500 = 0 * thetas;

  for k = 1 : 9
    theta = thetas(:, k);
    for mu = linspace((k * pi)^2 * 1.1, 5000, 100)
      [theta, history] = zeroOfBuck(length(theta), mu, 1e-6, theta);
    end
    at_mu_500(:, k) = theta;

    fst_theta = [0; theta; 0].';
    x = cumsum(sin(fst_theta)) / (n + 1);
    y = cumsum(cos(fst_theta)) / (n + 1);
    plot(fst_theta); hold on;
  end


end
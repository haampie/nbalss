function check_difference(k)
  A1 = dlmread('07_re_40_par19_0.0/fort.40', '', 2, 0);
  A2 = dlmread('10_re40_par_19_0.0_onderaan/fort.40', '', 2, 0);

  B1 = reshape(A1(:, k), 128, 128);
  B2 = reshape(A2(:, k), 128, 128);

  sum_of_solutions = flipud(B2) + B1;

  fprintf('Norm of flipud(B2) + B1 = %e\n', norm(sum_of_solutions, 'fro')); 

  surf(flipud(B2) + B1);
end
function plot_it(filename, k)
  A = dlmread(filename, '', 2, 0);
figure;
  contour(reshape(A(:, k), 128, 128));
end  
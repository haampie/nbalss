
  reynolds = [
   0.291449E+02  0.292888E+02 0.294318E+02 0.295739E+02 0.296528E+02 
  ];

  eigenvalues = [
      -0.1585E+00  +  -0.6932E-08i ;% -0.5470E+01  +  -0.9913E-07i -0.2156E+02  +   0.2367E-05i -0.2993E+02  +   0.1419E-03i;
       -0.9538E-01  +  -0.4786E-10i;% -0.5400E+01  +  -0.3081E-09i -0.2136E+02  +   0.1371E-06i -0.2977E+02  +   0.5216E-04i;
       -0.3375E-01  +  -0.1990E-10i ;%  -0.5331E+01  +  -0.1440E-08i   -0.2116E+02  +   0.6018E-07i  -0.2962E+02  +  -0.1850E-04i;
        0.2639E-01  +   0.9408E-09i;% -0.5262E+01  +  -0.1595E-07i -0.2096E+02  +  -0.5125E-05i -0.2946E+02  +   0.5740E-06i;
         0.5932E-01  +   0.2336E-06i ;% -0.5225E+01  +  -0.8771E-09i  -0.2086E+02  +   0.6621E-06i  -0.2938E+02  +   0.1273E-05i  
      ];

  figure;

  legend_labels = {};

  markers = {'+', '*', 'o', 'x'};

  marker_no = 1;

  for i = 1 : 1 : length(reynolds)
    scatter(real(eigenvalues(i, :)), imag(eigenvalues(i, :)), markers{marker_no}); hold on
    legend_labels{length(legend_labels) + 1} = [sprintf('Re = %2.2f', reynolds(i))];
    marker_no = marker_no + 1;
    if marker_no > length(markers)
      marker_no = 1;
    end
  end

  grid;

  ax = gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
  legend(legend_labels, 'Location','northwest')

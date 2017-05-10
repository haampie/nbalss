function plotQG(directory)

  figure;


  % resolution 
  n = 128;
  m = 128;
  % read solution
  [lab icp par xl xlp sig sol solup soleig]  = readfort3(n,m,[directory 'fort.3']);
  % constants 
  udim = 1.6e-02; 
  ldim = 1.0e+06; 
  hdim = 6.0e+02; 
  fact = udim*hdim*ldim/1.0e+06; 
  % grid
  for i=1:n
     x(i) = (i-1)/(n-1);
  end
  for j=1:m
     y(j) = (j-1)/(m-1); 
  end
  % scaling
  maxp = max(max(sol(:,:,2))) 
  % plot
  plot2D_zs(x,y,fact*sol(:,:,2)/maxp,10,1); 
  xlabel('x/L','FontSize',20)
  ylabel('y/L','FontSize',20)
  title('')
  saveeps(15,10, strcat(directory, 'Fig3d.eps'))
end
 




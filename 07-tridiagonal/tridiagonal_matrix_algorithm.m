function x = tridiagonal_matrix_algorithm(a, b, c, d)
  np = length(b);
  x  = zeros(1, np);
  
  c(1) = c(1)/b(1);
  d(1) = d(1)/b(1);
  
  for i = 2:np-1
    c(i) = c(i)/(b(i) - c(i-1) * a(i-1));
    d(i) = (d(i) - d(i-1) * a(i-1))/(b(i) - c(i-1) * a(i-1));
  end %for
  
  x(np) = (d(np) - d(np-1) * a(np-1))/(b(np) - c(np-1) * a(np-1));
  
  for i = np-1:-1:1
    x(i) = d(i) - c(i) * x(i+1);
  end %for
end
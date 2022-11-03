function [x, uan, unu] = getu(a, b, x0, xmax, np)
  h      = (xmax - x0) / (np - 1);
  x      = linspace(x0,xmax,np);
  uan    = a * cos(x) + b * sin(x);
  unu    = zeros(1, np); 
  unu(1) = a;
  unu(2) = h * b + a;

  for i=3:np
    unu(i) = (2 - h * h) * unu(i-1) + - unu(i-2);
  end
end
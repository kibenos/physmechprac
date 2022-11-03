function [x, uan, unu, unu_half] = getu(a, xmax, nh)
  x0     = 0;
  h      = (xmax - x0) / nh;
  x      = linspace(x0,xmax,nh+1);
  x_half = linspace(x0+h/2,xmax+h/2,nh+1);
  
  uan      = a * (2 - cos(x));
  unu      = zeros(1, nh + 1);
  unu_half = zeros(1, nh + 1);
  
  unu(1)      = a;
  unu_half(1) = a;
  
  for i=2:nh+1
    unu(i) = unu(i-1) + a * sin(x(i-1)) * h;
    unu_half(i) = unu_half(i-1) + a * sin(x_half(i-1)) * h;
  end
end
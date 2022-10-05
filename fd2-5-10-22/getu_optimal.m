function [x, uan, unu] = getu_optimal(need_numeric, order, a, b, x0, xmax, np)
  h   = (xmax - x0) / (np - 1);
  x   = linspace(x0,xmax,np);
  uan = a * cos(x) + b * sin(x);

  if need_numeric
      unu    = zeros(1, np); 
      unu(1) = a;
      unu(2) = a * cos(x0 + h) + b * sin(x0 + h);
    
      k   = 2:2:order;
      one = ones(1, length(k));
      sgn = (-one).^(k/2 - one);
      s   = h * one;
      s   = s.^k;
      s   = s./factorial(k);
      s   = sum(s.*sgn);

      c1 = 0.5 / s; % 1 / (h^2 - h^4 / 12);
      c0 = - 2 * c1;
      c  = - (1 + c0) / c1;
    
      for i=3:np
        unu(i) = c * unu(i-1) - unu(i-2);
      end
  else
      unu = 0;
  end
end
function [x, uan, unu] = getu_optimal(need_numeric, order, a, b, x0, xmax, np)
  h   = (xmax - x0) / (np - 1);
  x   = linspace(x0,xmax,np);
  uan = a * cos(x) + b * sin(x);

  if need_numeric    
      c1 = -0.5 / get_series_sum(false, order, h);
      c0 = - 2 * c1;
      c  = - (1 + c0) / c1;

      unu    = zeros(1, np); 
      unu(1) = a;
      unu(2) = unu(1) + a * get_series_sum(false, order, h) + b * get_series_sum(true, order, h);

      for i=3:np
        unu(i) = c * unu(i-1) - unu(i-2);
      end
  else
      unu = 0;
  end
end
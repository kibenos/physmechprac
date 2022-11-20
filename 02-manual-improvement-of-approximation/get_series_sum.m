function [s] = get_series_sum(is_odd, order, h)
    if is_odd
      k   = 1:2:order;
      one = ones(1, length(k));
      sgn = (-one).^((k+1)/2 - one);
      s   = h * one;
      s   = s.^k;
      s   = s./factorial(k);
      s   = sum(s.*sgn);
    else
      k   = 2:2:order;
      one = ones(1, length(k));
      sgn = (-one).^(k/2);
      s   = h * one;
      s   = s.^k;
      s   = s./factorial(k);
      s   = sum(s.*sgn);
    end
end

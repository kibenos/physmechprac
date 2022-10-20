clc
clf 
colormap jet
figure(1)

x0     = 0;
xmax   = 8 * pi;
diffop = [0 1 0];
itmax  = 10;
itbase = (sqrt(5)+1)/2;


%% init
x_prev     = [];
du_prev    = [];
eps_prev   = 0;
np_prev    = 0;
hmax_prev  = 0;

xfin  = linspace(x0, xmax, 2^(meshmax + 2));
dufin = danalytics(xfin);

%% table header
w  = max(7, round(log10(2^(meshmax + 2) + 1) + 1));
fprintf('%*cpoints |    max h     |      ε       |    rat ε    |    order    |\n', w - 6, '');
fprintf('%-*c-------|--------------|--------------|-------------|-------------|\n', w - 6, '-');

%% main cycle
for it=1:itmax
  nr    = ceil(itbase^(it + 1));
  ng    = ceil(itbase^(it + 1));

  %x     = linspace(x0, xmax, 2^(it+2)); 
  x     = makemesh(x0, xmax, x_prev, du_prev, nr, ng, 1.0e-17);
  np    = length(x);
  u     = analytics(x);
  du    = zeros(1, length(x));
  duan  = danalytics(x);

  for i=1:np
      if i == 1
          c     = getcoeff(diffop, [0 x(2)-x(1) x(3)-x(1)]);
          du(i) = u(1:3) * c;
      elseif i == length(x)
          c     = getcoeff(diffop, [x(end-2)-x(end) x(end-1)-x(end) 0]);
          du(i) = u(end-2:end) * c;
      else
          c     = getcoeff(diffop, [x(i-1)-x(i) 0 x(i+1)-x(i)]);
          du(i) = u(i-1:i+1) * c;
      end
  end

  eps   = max(abs(duan - du)); % sqrt(trapz(x, (duan - du).^2));
  hmax  = max(x(2:end) - x(1:end-1));
  outlg = sprintf('%d points; max(h) = %e; ε = %e;', np, hmax, eps);
  outcl = sprintf('%*d | %e | %e |', w, np, hmax, eps);

  if it > 1
    rateps = eps_prev/eps;
    order  = log(rateps) / log(hmax_prev/hmax);
    outlg  = strcat(outlg, sprintf(' rat ε = %f; order = %f; ', rateps, order));
    outcl  = strcat(outcl, sprintf(' %11.6f | %11.6f |', rateps, order));
  else
      outcl = strcat(outcl, '             |             |');
  end
  
%% plot&out
  clf;
  hold on;
  xlabel('x');
  ylabel('y');
  plot(xfin, dufin, 'color', 'green');
  if it > 1
    plot(x_prev, du_prev, '--.', 'color', [.6, .6, .6]);
    plot(x, du, '--.', 'color',' red');
    legend('analytics', sprintf('%d points', np_prev), sprintf('%d points', np));
  else
    plot(x, du, '--.', 'color', 'red');
    legend('analytics', sprintf('%d points', np  ));
  end
  title(outlg);
  
  outcl = strcat(outcl, '\n') ;
  fprintf(outcl);
  pause;

%% save prev
  x_prev    = x;
  du_prev   = du;
  eps_prev  = eps;
  np_prev   = np;
  hmax_prev = hmax;
end








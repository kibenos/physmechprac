clc
clf 
colormap jet
figure(1)

x0     = 1;
xmax   = 8 * pi;
npinit = 32;
itmax  = 10;
meps   = 1e-4;

diffop = [ 0  0  1];
pat    = []; % [-1  0  1];
lpat   = []; % [ 0  1  2];
rpat   = []; % [-2 -1  0];

alpha       = 0.95;
adapt_order = 3;    % order of the derivative influencing on adaptive meshing 

%% init
x_prev    = [];
u_prev    = [];
du_prev   = [];
duan_prev = [];
eps_prev  = 0;
np_prev   = 0;
hmax_prev = 0;
hmin_prev = 0;

xfin  = linspace(x0, xmax, 100);
dufin = danalytics(xfin);

%% table header
w  = max(7, round(log10(2^(itmax + 2) + 1) + 1));
fprintf('%spoints |    min h     |    max h     |     eps      |   rat eps   |\n', repmat(' ', 1, w - 6));
fprintf('%s-------|--------------|--------------|--------------|-------------|\n', repmat('-', 1, w - 6));

%% main cycle
for it=1:itmax
  
% init random with adaptive meshing
  if it == 1
      x    = linspace(x0, xmax, npinit);
      u    = analytics(x);
      duan = danalytics(x);
  else
      diffop_adapt      = zeros(1, adapt_order + 1);
      diffop_adapt(end) = 1;
      du_adapt          = fdiff(x_prev, u, diffop_adapt, [], [], []);
      du_adapt_norm     = sqrt(sum(du_adapt.^2))/length(du_adapt);
      wraps             = zeros(1, np_prev);
     
      for i = 1:np_prev
          if i == 1
              h = x_prev(2) - x_prev(1);
          elseif i == np_prev
              h = x_prev(end) - x_prev(end-1);
          else
              h = min(x_prev(i+1) - x_prev(i), x_prev(i) - x_prev(i-1));
          end

          wraps(i) = h / (alpha * abs(du_adapt(i)) + (1 - alpha) * du_adapt_norm);
      end

      [x, u, duan] = wrapmesh(x_prev, u_prev, duan_prev, wraps, pat, lpat, rpat, length(diffop)-1, meps);
      
      % fix this shit
      u    = analytics(x);
      duan = danalytics(x);
  end

  np = length(x);
  du = fdiff(x, u, diffop, pat, lpat, rpat);

%   eps   = sqrt(trapz((duan - du).^2));
%   eps   = sqrt(sum((duan - du).^2))/length(du); % L2(x0, xmax) norm
  eps   = max(abs(duan - du));                  %  C[x0, xmax] norm
  h     = x(2:end) - x(1:end-1);
  hmax  = max(h);
  hmin  = min(h);

  outlg = sprintf('%d points; min(h) = %e; max(h) = %e; eps = %e;', np, hmin, hmax, eps);
  outcl = sprintf('%*d | %e | %e | %e |', w, np, hmin, hmax, eps);

  if it > 1
    rateps = eps_prev/eps;
    outlg  = strcat(outlg, sprintf(' rat eps = %f; ', rateps));
    outcl  = strcat(outcl, sprintf(' %11.6f |', rateps));
  else
      outcl = strcat(outcl, '             |');
  end
  
%% plot&out
  clf;
  hold on;
  xlabel('x');
  ylabel('y');
  plot(xfin, dufin, 'color', 'green');
  if it > 1
    plot(x_prev, du_prev, '--.', 'color', [.6, .6, .6]);
    plot(x, du, '--.', 'color', 'red');
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
  u_prev    = u;
  du_prev   = du;
  duan_prev = duan;
  eps_prev  = eps;
  np_prev   = np;
  hmax_prev = hmax;
  hmin_prev = hmin;
end








clc
clf 
colormap jet
figure(1)

x0     = 0;
xmax   = 20;
diffop = [0 1];
itmax  = 10;
itbase = 2;

%% seting pattern
patsz     = 3; % pattern size
cpattern  = zeros(1, patsz);
lpatterns = zeros(floor(patsz/2), patsz);
rpatterns = zeros(floor(patsz/2), patsz);

if mod(patsz, 2)
  cpattern = ceil(patsz/2)-patsz:floor(patsz/2)
else
  cpattern = [-patsz/2:-1, 1:patsz/2]
end

for i=1:floor(patsz/2)
  lpatterns(i, :) = 1-i:patsz-i;
  rpatterns(i, :) = -patsz+i:i-1;
end

lpatterns
rpatterns = rpatterns(end:-1:1,:)

%% init
x_prev    = [];
du_prev   = [];
eps_prev  = 0;
np_prev   = 0;
hmax_prev = 0;
hmin_prev = 0;

xfin = linspace(x0, xmax, 200);
duanfin = danalytics(xfin);

%% table header
w  = max(7, round(log10(itbase^(itmax + 2) + 1) + 1));
fprintf('%spoints |    min h     |    max h     |     eps      |   rat eps   |    order    |\n', repmat(' ', 1, w - 6));
fprintf('%s-------|--------------|--------------|--------------|-------------|-------------|\n', repmat('-', 1, w - 6));

%% main cycle
for it=1:itmax
  x    = linspace(x0, xmax, ceil(itbase^(it + 2)));
  np   = length(x);
  u    = analytics(x);
  du   = zeros(1, length(x));
  duan = danalytics(x);
  
  D  = getmatrix(diffop, x, cpattern, lpatterns, rpatterns);
  du = u*D';

%   eps   = sqrt(trapz((duan - du).^2));
%   eps   = sqrt(sum((duan - du).^2))/length(du); % L2(x0, xmax) norm
  eps   = max(abs(duan - du));                  %  C[x0, xmax] norm
  h     = x(2:end) - x(1:end-1);
  hmax  = max(h);
  hmin  = min(h);

  outlg = sprintf('%d points; max(h) = %e; min(h) = %e; eps = %e;', np, hmin, hmax, eps);
  outcl = sprintf('%*d | %e | %e | %e |', w, np, hmin, hmax, eps);

  if it > 1
    rateps = eps_prev/eps;
    order  = log(rateps) / log(hmax_prev/hmax);
    outlg  = strcat(outlg, sprintf(' rat eps = %f; order = %f; ', rateps, order));
    outcl  = strcat(outcl, sprintf(' %11.6f | %11.6f |', rateps, order));
  else
      outcl = strcat(outcl, '             |             |');
  end
  
%% plot&out
  clf;
  hold on;
  xlabel('x');
  ylabel('y');
  plot(xfin, duanfin, '--.', 'color', 'green');
  if it > 1
    plot(x_prev, du_prev, '--.', 'color', [.6, .6, .6]);
    plot(x, du, '--.', 'color', 'red');
    legend('analytics', sprintf('%d points', np_prev), sprintf('%d points', np));
  else
    plot(x, du, '--.', 'color', 'red');
    legend('analytics', sprintf('%d points', np));
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
  hmin_prev = hmin;
end








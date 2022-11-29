clc
clf 
colormap jet
figure(1)

x0     = 0;
xmax   = 2 * pi;
diffop = [1 0 1];

a      = 1;
b      = 0;

itmax  = 10;
itbase = 2;

%% symmetric pattern
%schord    = 2;
%patsz     = schord - 1 + find(diffop, 1, 'last'); % pattern size
%cpattern  = zeros(1, patsz);
%lpatterns = zeros(floor(patsz/2), patsz);
%rpatterns = zeros(floor(patsz/2), patsz);
%
%if mod(patsz, 2)
%  cpattern = ceil(patsz/2)-patsz:floor(patsz/2);
%else
%  cpattern = [-patsz/2:-1, 1:patsz/2];
%end
%
%for i=1:floor(patsz/2)
%  lpatterns(i, :) = 1-i:patsz-i;
%  rpatterns(i, :) = -patsz+i:i-1;
%end
%
%rpatterns = rpatterns(end:-1:1,:);

%% pattern
cpattern  = [-1  0  1];
lpatterns = [ 0  1  2];
rpatterns = [-2 -1  0];

%% init
x_prev    = [];
u_prev    = [];
eps_prev  = 0;
np_prev   = 0;
hmax_prev = 0;
hmin_prev = 0;

%% plotting analytics
xfin   = linspace(x0, xmax, 200);
uanfin = analytics(xfin, a, b);

%% table header
w  = max(7, round(log10(itbase^(itmax + 2) + 1) + 1));
fprintf('%spoints |    min h     |    max h     |     eps      |   rat eps   |    order    |\n', repmat(' ', 1, w - 6));
fprintf('%s-------|--------------|--------------|--------------|-------------|-------------|\n', repmat('-', 1, w - 6));

%% main cycle
for it=1:itmax
  x   = linspace(x0, xmax, ceil(itbase^(it + 2)));
  np  = length(x);
  uan = analytics(x, a, b);
  D   = getmatrix(diffop, x, cpattern, lpatterns, rpatterns);

  %% RHS
  RHS = rhs(x)';
  
  %% BC
  BCdifop       = [0 1];
  BCpattern     = [-1 0];
  
  D(3:end, :)   = D(2:end-1, :);
  D(1:2, :)     = zeros(2, np);
  D(1, 1)       = 1;
  D(2, 1:2)     = getcoeff(BCdifop, x(2 + BCpattern) - x(2));
  
  RHS(3:end, :) = RHS(2:end-1, :);
  RHS(1:2, :)   = zeros(2, 1);
  RHS(1, 1)     = a;
  RHS(2, 1)     = b;
  
  %% SOLUTION
  u = linsolve(D, RHS)';

  %% ERROR
  %eps   = sqrt(sum((uan - u).^2))/length(u);  % L2(x0, xmax) norm
  eps   = max(abs(uan - u));                  %  C[x0, xmax] norm
  h     = x(2:end) - x(1:end-1);
  hmax  = max(h);
  hmin  = min(h);

  %% output formatting
  outleg = sprintf('%d points; max(h) = %e; min(h) = %e; eps = %e;', np, hmin, hmax, eps);
  outcon = sprintf('%*d | %e | %e | %e |', w, np, hmin, hmax, eps);

  if it > 1
    rateps = eps_prev/eps;
    order  = log(rateps) / log(hmax_prev/hmax);
    outleg  = strcat(outleg, sprintf(' rat eps = %f; order = %f; ', rateps, order));
    outcon  = strcat(outcon, sprintf(' %11.6f | %11.6f |', rateps, order));
  else
      outcon = strcat(outcon, '             |             |');
  end
  
  %% plotting
  clf;
  hold on;
  xlabel('x');
  ylabel('y');
  plot(xfin, uanfin, '--.', 'color', 'green');
  if it > 1
    plot(x_prev, u_prev, '--.', 'color', [.6, .6, .6]);
    plot(x, u, '--.', 'color', 'red');
    legend('analytics', sprintf('%d points', np_prev), sprintf('%d points', np));
  else
    plot(x, u, '--.', 'color', 'red');
    legend('analytics', sprintf('%d points', np));
  end
  title(outleg);
  
  %% output
  outcon = strcat(outcon, '\n') ;
  fprintf(outcon);

  %% save prev
  x_prev    = x;
  u_prev    = u;
  eps_prev  = eps;
  np_prev   = np;
  hmax_prev = hmax;
  hmin_prev = hmin;
  
  if it < itmax
    pause;
  end
end








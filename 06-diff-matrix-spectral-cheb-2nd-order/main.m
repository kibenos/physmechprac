clc
clf 
colormap jet
figure(1)

x0     = -1;
xmax   = 1;
diffop = [1 0 1];

a      = 2;
b      = 3;

itmax  = 27;
itbase = 2;

%% init
x_prev    = [];
u_prev    = [];
eps_prev  = 0;
np_prev   = 0;
hmax_prev = 0;
hmin_prev = 0;

hold on;
xlabel('N');
ylabel('error');

%% table header
w  = max(7, round(log10(itbase^(itmax + 2) + 1) + 1));
fprintf('%spoints |     min h     |    max h      |     eps      |   rat eps   |    order    |\n', repmat(' ', 1, w - 6));
fprintf('%s-------|---------------|---------------|--------------|-------------|-------------|\n', repmat('-', 1, w - 6));

%% main cycle
for it=1:itmax
  [D1, x] = getmatrix_cheb(it + 2);
  np      = length(x);
  D2      = D1 * D1;
  D0      = eye(np);
  uan     = analytics(x, a, b);
  
  %% SYSTEM
  M = D2 + D0;

  %% RHS
  RHS = rhs(x);
  
  %% BC
  M(1, :)     = zeros(1, np);
  M(1, 1)     = 1;
  M(end, :)   = D1(end, :);
  RHS(1, 1)   = analytics(x(1), a, b);
  RHS(end, 1) = danalytics(x(end), a, b);
  
  %% SOLUTION
  u = linsolve(M, RHS);

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
    rateps  = eps_prev/eps;
    order   = log(rateps) / log(hmax_prev/hmax);
    outleg  = strcat(outleg, sprintf(' rat eps = %f; order = %f; ', rateps, order));
    outcon  = strcat(outcon, sprintf(' %11.6f | %11.6f |', rateps, order));
  else
      outcon = strcat(outcon, '             |             |');
  end
  
  %% plotting
  plot(log10(np), log10(eps), '.', 'color', 'red');
  
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








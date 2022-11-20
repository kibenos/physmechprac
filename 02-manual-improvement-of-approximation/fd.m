clc
clf 
colormap jet

a       = 0;
b       = 1;
x0      = 0;
xmax    = 2 * pi;
meshmax = 9;
order   = 4;

%% init
x          = [];
x_prev     = [];
u          = []; 
u_prev     = [];
uan        = []; 
eps_prev   = 0;
delta_prev = 0;
np_prev    = 0; 
[xfin, uanfin, ~] = getu_optimal(false, order, a, b, x0, xmax, 2^(meshmax + 2));

%% table header
w  = max(7, round(log10(2^(meshmax + 2) + 1) + 1));
wo = round(log10(order^2)) + 1;
fprintf('%*cpoints |      h       |      ε       |      δ       |    rat ε    |    rat δ\n', w - 6, '');
fprintf('%-*c-------|--------------|--------------|--------------|-------------|------------\n', w - 6, '-');

%% main cycle
for mesh=1:meshmax
  np          = 2^(mesh + 2) + 1;
  [x, uan, u] = getu_optimal(true, order, a, b, x0, xmax, np);
  eps         = max(abs(uan - u)); % sqrt(trapz(x, (uan - u).^2));

  outlg = sprintf('%dth order accurate scheme; %d points; h = %e; ε = %e;', order, np, (xmax - x0) / (np - 1), eps);
  outcl = sprintf('%*d | %e | %e |', w, np, (xmax - x0) / (np - 1), eps);
  
  if mesh > 1
    delta  = max(abs(u_prev - u(1:2:end))); % sqrt(trapz(x_prev, (u_prev - u(1:2:end)).^2));
    rateps = eps_prev/eps;

    outlg = strcat(outlg, sprintf(' δ = %e; rat ε = %f;', delta, rateps));
    outcl = strcat(outcl, sprintf(' %e | %11.6f |', delta, rateps));
    if mesh > 2
        ratdelta = delta_prev/delta;
        outlg = strcat(outlg, sprintf(' rat δ = %f', ratdelta));
        outcl = strcat(outcl, sprintf(' %11.6f', ratdelta));
    end

    delta_prev = delta;
  else
      outcl = strcat(outcl, '              |             |');
  end
  
%% plot&out
  clf;
  hold on;
  xlabel('x');
  ylabel('y');
  plot(xfin, uanfin, 'color', 'green');
  if mesh > 1
    plot(x_prev, u_prev, '--.', 'color', [.6, .6, .6]);
    plot(x, u, '--.', 'color',' red');
    legend('analytics', sprintf('%d points', np_prev), sprintf('%d points', np));
  else
    plot(x, u, '--.', 'color', 'red');
    legend('analytics', sprintf('%d points', np  ));
  end
  title(outlg);
  
  outcl = strcat(outcl, '\n') ;
  fprintf(outcl);
  pause;

%% save prev
  u_prev   = u;
  x_prev   = x;
  eps_prev = eps;
  np_prev  = np;
end








a       = 1;
xmax    = 5;
meshmax = 5;

x           = [];
x_prev      = [];
uan         = []; 
u           = []; 
u_half      = [];
u_prev      = [];
eps_prev    = 0;
delta_prev  = 0;

is_half_scheme = false;

xlabel("x");
ylabel("y");

[xfin, uanfin, _, _] = getu(a, xmax, 2^(meshmax + 2));

for mesh=1:meshmax
  nh = 2^(mesh + 2);
  
  if is_half_scheme
    [x, uan, _, u] = getu(a, xmax, nh);
  else
    [x, uan, u, _] = getu(a, xmax, nh);
  end
  
  eps = max(abs(uan - u));
  out = sprintf("M = %d; \\epsilon = %f; ", nh + 1, eps);
  
  if mesh > 1
    delta = max(abs(u_prev - u(1:2:end)));
    out = [out, sprintf("\\delta = %f; rat \\epsilon = %f; ", delta, eps_prev/eps)];
    
    if mesh > 2
      out = [out, sprintf("rat \\delta = %f", delta_prev/delta)];
    end
    
    delta_prev = delta;
  end
  
  clf;
  hold on;
  plot(xfin,uanfin,"color","green");
  plot(x,u,"--.");
  if mesh > 1
    plot(x_prev,u_prev,"--.");
    legend("analytics", sprintf("M = %d", nh), sprintf("M = %d", 2^(mesh + 1)));
  else
    legend("analytics", sprintf("M = %d", nh));
  end
  title(out);
  
  printf([out, "\n"]);
  pause;
  
  u_prev = u;
  x_prev = x;
  eps_prev = eps;
end








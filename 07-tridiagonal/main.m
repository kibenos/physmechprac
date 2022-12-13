clc
clf 
colormap jet
%figure(1)

itmax  = 9;
itbase = 2;
meps   = 1e-1;

%hold on;
%xlabel('N');
%ylabel('error');

%% table header
w  = max(7, round(log10(itbase^(itmax + 2) + 1) + 1));
fprintf('%spoints |     eps      |\n', repmat(' ', 1, w - 6));
fprintf('%s-------|--------------|\n', repmat('-', 1, w - 6));

%% main cycle
for it=1:itmax
  np = ceil(itbase^(it + 2));
  a  = -1 * ones(1, np-1) + meps * rand(1, np-1);
  b  =  3 * ones(1, np  ) + meps * rand(1, np  );
  c  = -1 * ones(1, np-1) + meps * rand(1, np-1);
  d  = zeros(1, np);
  M  = diag(a, -1) + diag(b) + diag(c, 1);
  
  x1 = tridiagonal_matrix_algorithm(a, b, c, d);
  x2 = linsolve(M, d')';
  
  eps = max(abs(x1 - x2));

  %% output formatting
  outleg = sprintf('%d points; eps = %e;', np, eps);
  outcon = sprintf('%*d | %e |', w, np, eps);
  
  %% plotting
  %plot(log10(np), log10(eps), '.', 'color', 'red');
  
  %% output
  outcon = strcat(outcon, '\n') ;
  fprintf(outcon);
  
  if it < itmax
    pause;
  end
end








% diffop    - array of differential operator coefficients
%        example:
%             diffop:
%             [1 0 -2]
%             operator:
%             (I - 2 * d^2)[u(x)] = u(x) - 2 * (d^2/dx^2)u(x)
%
% mesh      - array of poins
%        denote:
%        mesh(1:len(lpatterns))                    - left points
%        mesh(len(lpatterns)+1:end-len(rpatterns)) - inner points
%        mesh(end-len(rpatterns)+1:end)            - right points
%
% cpattern  - array of pattern indexes offsets for each inner point
%             (offsets from the point where the diffop is calculated)
%        example:
%             cpattern:
%             [-2 0 1]
%             pattern: 
%             ---*---.---*---*---
%                        ^ derivative at this point
%
% lpatterns - array of patterns for left points
%        pattern lpatterns(i, :) correspond mesh(i) point
%
% rpatterns - array of patterns for right points
%        pattern rpatterns(i, :) correspond mesh(end-len(rpatterns)+i) point

function D = getmatrix(diffop, mesh, cpattern, lpatterns, rpatterns)
    np = length(mesh);
    D  = sparse(np, np);
    
    for i = 1:np
      % left points
      if i <= size(lpatterns, 1)
          D(i, i+lpatterns(i, :)) = getcoeff(diffop, mesh(i+lpatterns(i, :)) - mesh(i));
      
      % inner points
      elseif i <= np-size(rpatterns, 1)
          D(i, i+cpattern) = getcoeff(diffop, mesh(i+cpattern) - mesh(i));
      
      % right points
      elseif i > np-size(rpatterns, 1) 
          D(i, i+rpatterns(end+i-np, :)) = getcoeff(diffop, mesh(i+rpatterns(end+i-np, :)) - mesh(i));
      end
    end
end
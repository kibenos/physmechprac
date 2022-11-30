% diffop  - array of differential operator coefficients
%        example:
%             diffop:
%             [1 0 -2]
%             operator:
%             (I - 2 * d^2)[u(x)] = u(x) - 2 * (d^2/dx^2)u(x)
%
% offsets - array of pattern offsets for each inner point
%             (derivative is calculated at a point with zero offset)
%        example:
%             offsets:
%             [-1.0 0.0 0.5]
%             pattern: 
%             ----*-------.--------*-------*-----> 
%             ---(-1)---(-0.5)---(0.0)---(0.5)---> x
%                                  ^ derivative at this point

function c = getcoeff(diffop, offsets)

    if length(diffop) < length(offsets)
        diffop = cat(2, diffop, zeros(1, length(offsets) - length(diffop)));
    elseif length(diffop) > length(offsets)
        fprintf("error: diff operator exceeds scheme size\n");
        return;
    end

    A = zeros(length(diffop), length(offsets));

    for i=1:length(diffop)
        for j=1:length(offsets)
            A(i, j) = (offsets(j))^(i-1)/factorial(i-1);
        end
    end

    c = linsolve(A, diffop');
end
function [c, A] = getcoeff(diff_coeffs, sch_steps)

    A = zeros(length(diff_coeffs), length(sch_steps));

    for i=1:length(diff_coeffs)
        for j=1:length(sch_steps)
            A(i, j) = (sch_steps(j))^(i-1)/factorial(i-1);
        end
    end

    c = linsolve(A, diff_coeffs');
end
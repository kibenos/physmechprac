function c = getcoeff(diff_coeffs, sch_steps)

    if length(diff_coeffs) < length(sch_steps)
        diff_coeffs = cat(2, diff_coeffs, zeros(1, length(sch_steps) - length(diff_coeffs)));
    elseif length(diff_coeffs) > length(sch_steps)
        fprintf("error: diff operator exceeds scheme size\n");
        return;
    end

    A = zeros(length(diff_coeffs), length(sch_steps));

    for i=1:length(diff_coeffs)
        for j=1:length(sch_steps)
            A(i, j) = (sch_steps(j))^(i-1)/factorial(i-1);
        end
    end

    c = linsolve(A, diff_coeffs');
end
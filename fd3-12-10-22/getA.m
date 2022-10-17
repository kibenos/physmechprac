function [A] = getA(diff_coeff, sch_ordinal, h)
    A = zeros(length(diff_coeff), length(sch_ordinal));
    for i=1:length(diff_coeff)
        for j=1:length(sch_ordinal)
            A(i, j) = (h * sch_ordinal(j))^(i-1)/factorial(i-1);
        end
    end
end
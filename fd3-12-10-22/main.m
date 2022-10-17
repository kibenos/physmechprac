%format rational

diff_coeff  = [ 0  0 1];
sch_ordinal = [-2 -1 0];
h           = 1/10;

A = getA(diff_coeff, sch_ordinal, h);
c = linsolve(A, diff_coeff');

text = [];
for i=1:length(sch_ordinal)
    text = [text num2str(c(i, 1)) ' \cdot u_{' int2str(sch_ordinal(i)) '} '];
    if i < length(sch_ordinal)
        text = [text '+ '];
    end
end

title(text);
function [x, u, duan] = splitmesh(x_prev, u_prev, duan_prev, factors)
    split = {[]};
    idx   = {[]};
    
    for i = 1:length(x_prev)-1
        split{i}  = linspace(x_prev(i), x_prev(i+1), factors(i)+1);
        split{i}  = split{i}(1:end-1);

        idx{i}    = zeros(1, factors(i));
        idx{i}(1) = 1;
    end

    x          = cat(2, split{:});
    x(end+1)   = x_prev(end);
    idx        = logical(cat(2, idx{:}));
    idx(end+1) = 1;
    
    u          = zeros(1, length(x));
    duan       = zeros(1, length(x));

    u(idx)     = u_prev;
    duan(idx)  = duan_prev;

    u(~idx)    = analytics(x(~idx));
    duan(~idx) = danalytics(x(~idx));
end
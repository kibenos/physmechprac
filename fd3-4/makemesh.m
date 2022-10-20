function mesh = makemesh(a, b, auxmesh, auxval, nr, ng, eps)
    x  = [];
    if isempty(auxval) && nr ~= 0
        x = cat(2, x, (b - a).*rand(1, nr) + a);
    elseif isempty(auxval) && ng ~= 0
        x = cat(2, x, normrnd((a+b)/2, (b-a)/6, 1, ng));
    elseif ~isempty(auxval)
        x0 = a;
        m  = median(abs(auxval));

        is_uniform = auxval(1) < m;
        for i=2:length(auxval)
            if is_uniform && (abs(auxval(i)) > m || i == length(auxval))
                x = cat(2, x, makemesh(x0, auxmesh(i), [], [], nr, 0, eps));
                x0 = auxmesh(i);
                is_uniform = false;
            elseif ~is_uniform && (abs(auxval(i)) < m || i == length(auxval))
                x = cat(2, x, makemesh(x0, auxmesh(i), [], [], 0, ng, eps));
                x0 = auxmesh(i);
                is_uniform = true;
            end
        end
    end

    x = sort(x);
    mesh = zeros(1, length(x) + 2);
    mesh(1) = a;
    
    midx = 1;
    for i=1:length(x)
        if x(i) - mesh(midx) > eps && b - x(i) > eps
            midx = midx + 1;
            mesh(midx) = x(i);
        end
    end

    mesh(midx+1) = b;
    mesh = mesh(1:midx+1);
end
function [x, u, duan] = wrapmesh(x_prev, u_prev, duan_prev, wraps, cpattern, lpatterns, rpatterns, order, eps)
    wrap = {[]};
    idx  = {[]};

    [cpattern, lpatterns, rpatterns] = getdefaultpatterns(cpattern, lpatterns, rpatterns, order);

    np = length(x_prev);
    for i = 1:np
        if i <= size(lpatterns, 1)
            pattern = lpatterns(i,:);
        elseif i == size(lpatterns, 1)+1
            pattern = cpattern;
        elseif i > np-size(rpatterns, 1)
            pattern = rpatterns(np+1-i,:);
        end

        if wraps(i) > 0
            wrap{i}          = x_prev(i) + wraps(i).*pattern;
            idx{i}           = zeros(1, length(pattern));
            idx{i}(~pattern) = i;
        else
            wrap{i} = x_prev(i);
            idx{i}  = i;
        end
    end

    x   = cat(2, wrap{:});
    idx = cat(2, idx{:});

    % validate mesh
    mesh_prev  = [x; idx];
    mesh_prev  = sortrows(mesh_prev.', 1).'; % sort by first row
    mesh       = zeros(2, size(mesh_prev, 2) + 2);
    mesh(:, 1) = [x_prev(1); idx(1)];

    midx       = 1;
    for i=1:size(mesh_prev, 2)
        if mesh_prev(1, i) - mesh(1, midx) > eps && x_prev(end) - mesh_prev(1, i) > eps
            midx          = midx + 1;
            mesh(:, midx) = mesh_prev(:, i);
        elseif x_prev(end) - mesh_prev(1, i) > eps
            mesh(2, midx) = max(mesh(2, midx), mesh_prev(2, i));
        end
    end

    mesh(:, midx+1) = [x_prev(end); idx(end)];
    mesh            = mesh(:, 1:midx+1);

    x    = mesh(1, :);
    idx  = mesh(2, :);

    u    = zeros(1, length(x));
    duan = zeros(1, length(x));

    u(logical(idx))     = u_prev(nonzeros(idx)');
    duan(logical(idx))  = duan_prev(nonzeros(idx)');

    u(~idx)             = analytics(x(~idx));
    duan(~idx)          = danalytics(x(~idx));
end
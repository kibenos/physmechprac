function du = fdiff(mesh, u, diffop, cpattern, lpatterns, rpatterns)
    [cpattern, lpatterns, rpatterns] = getdefaultpatterns(cpattern, lpatterns, rpatterns, length(diffop)-1);
    
    np = length(mesh);
    du = zeros(1, np);

    for i = 1:np
        if i <= size(lpatterns, 1)
            pattern = lpatterns(i,:);
        elseif i == size(lpatterns, 1) + 1
            pattern = cpattern;
        elseif i > np-size(rpatterns, 1)
            pattern = rpatterns(np+1-i,:);
        end

        c     = getcoeff(diffop, mesh(pattern + i) - mesh(i));
        du(i) = u(pattern + i) * c;
    end

end
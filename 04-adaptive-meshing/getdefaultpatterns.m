function [cpattern, lpatterns, rpatterns] = getdefaultpatterns(cpat, lpats, rpats, order)

    nd  = length(cpat);  % number of points in inner pattern
    rnd = sum(cpat > 0); % number of right points in pattern
    lnd = sum(cpat < 0); % number of left  points in pattern

    % symmetrical default pattern at internal points
    if isempty(cpat)    
        rnd = ceil(order/2);
        lnd = rnd;
        nd  = 2 * rnd + 1;  

        cpattern = -rnd:lnd;
    else
        cpattern = cpat;
    end

    % default patterns near left boundary
    if isempty(lpats)
        lpatterns = zeros(lnd, nd);
        for i = 1:lnd
            lpatterns(i,:) = 1-i:nd-i;
        end
    else
        lpatterns = lpats;
    end

    % default patterns near right boundary
    if isempty(rpats)
        rpatterns = zeros(rnd, nd);
        for i = 1:rnd
            rpatterns(i,:) = -nd+i:i-1;
        end
    else
        rpatterns = rpats;
    end
end
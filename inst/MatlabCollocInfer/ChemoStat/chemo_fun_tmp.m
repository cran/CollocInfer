function r = chemo_fun_tmp(times, y, p, more)

    if nargin < 4, more = []; end

    r = y;
    p = exp(p);
    Q = p(3) .* y(:, 2) + p(4) .* y(:, 3);
    Qs = Q .* exp(10 * (Q - p(16)))./(1 + exp(10 .* (Q - p(16)))) + p(16)./(1 + exp(10 .* (Q - p(16))));

    r(:, 1) = p(6) .* (p(5) - y(:, 1)) - p(12) * y(:, 2) .* y(:, 1)./(p(10) + y(:, 1)) - ...
        p(12) .* y(:, 3) .* y(:, 1)./(p(11) + y(:, 1));
    
    r(:, 2) = y(:, 2) .* (p(9) .* p(12) .* y(:, 1)./(p(10) + y(:, 1)) - ...
        p(3) .* p(13) .* (y(:, 4) + y(:, 5))./(p(15) + Qs) - p(6));
    
    r(:, 3) = y(:, 3) .* (p(9) * p(12) .* y(:, 1)./(p(11) + y(:, 1)) - ...
        p(4) * p(13) .* (y(:, 4) + y(:, 5))./(p(15) +  Qs) - p(6));
    
    r(:, 4) = y(:, 4) .* (p(14) * p(13) * Q./(p(15) + Qs) - (p(6) + p(7) + p(8)));
    
    r(:, 5) = p(8) .* y(:, 4) - (p(6) + p(7)) .* y(:, 5);        
end

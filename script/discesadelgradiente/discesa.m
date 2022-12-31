
function [xk,lista_punti,kterm] = discesa(A, b, x0, nmax, toll,lista_punti)
    xk = x0;
    rk = b-A*xk;
    k = 0;
    while ((norm(rk) >= toll) && (k < nmax))
        lista_punti{k+1,1} = xk;
        dk = rk;
        alphak = dk' * rk / (dk' * A * dk);
        xk = xk + alphak * dk;
        rk = b - A * xk;
        k = k+1;
    end
    kterm = k;
    if k == nmax
        warning('Convergenza non raggiunta in nmax iterazioni');
    end
end


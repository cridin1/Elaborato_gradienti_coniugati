function [xk,lista_punti,kterm] = gradienteconiugato(A, b, x0, nmax, toll,lista_punti)
    xk = x0;
    rk = b-A*xk;
    vk = rk;
    k = 0;
    while ((norm(rk) >= toll) && (k < nmax))
        lista_punti{k+1,1} = xk;
        alphak = vk' * rk / (vk' * A * vk);
        xk = xk + alphak * vk;
        rk = b - A * xk;
        betak = vk' * A * rk / (vk' * A * vk);
        vk = rk-betak*vk;
        k = k+1;
    end
    lista_punti{k+1,1} = xk;
    kterm = k+1;
    if k == nmax
        warning('Convergenza non raggiunta in nmax iterazioni');
    end
end
function [xk,lista_punti,kterm] = p_gradienteconiugato(A, b, x0, nmax, toll, R1, R2, lista_punti)
    R = R1*R2;
    xk = x0;
    rk = b-A*xk;
    zk = R\rk;
    vk = zk;
    k = 0;
    while ((norm(rk) >= toll) && (k < nmax))
        lista_punti{k+1,1} = xk;
        alphak = vk' * rk / (vk' * A * vk);
        xk = xk + alphak * vk;
        rk = b - A * xk;
        zk = R\rk;
        betak = vk' * A * zk / (vk' * A * vk);
        vk = zk-betak*vk;
        k = k+1;
    end
    lista_punti{k+1,1} = xk;
    kterm = k+1;
    if k == nmax
        warning('Convergenza non raggiunta in nmax iterazioni');
    end
end
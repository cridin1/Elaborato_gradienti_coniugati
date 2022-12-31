
import discesa.*

N = 2;
mu = 6;
A = full(sprandsym(N, 1, 1/mu, 1));  %costruisco la matrice simmetrica e definita positiva 
                                      %(dim, densità, 1/indice_condizionamento, definita positiva = 1) 


%parametri
b = rand(N,1);
x0 = rand(N,1);       %considero il vettore nullo come pos iniziale
nmax = 50;
toll = 10^(-3);

%variabili per memorizzare dati
kterm = 0;
lista_punti = cell(nmax,1);


%algoritmo
[xk,lista_punti,kterm] = discesa(A, b, x0, nmax, toll,lista_punti);

%valutazione errore rispetto al metodo diretto
xt = A\b;
er = abs(xk-xt);
rk = b-A*xk;


%visualizzazione
syms x y
f = symfun(0.5*(A(1,1) * x^2 + (A(1,2)+A(2,1)) * x*y + A(2,2) * y^2) - b(1)*x - b(2)*y,[x y]);

%tuning dei parametri per la visualizzazione
camp = 0.01; %passo di campionamento
off = norm(x0); %scelgo offset
[X,Y] = meshgrid(-off+xt(1):camp:off+xt(1),-off+xt(2):camp:off+xt(2)); %griglia di punti


Z = 0.5*(A(1,1) * X.*X + (A(1,2)+A(2,1)) * X.*Y + A(2,2) * Y.*Y) - b(1)*X - b(2)*Y;
figure(1);
hold on


v = lista_punti{1,1};
xp = v(1);
yp = v(2);
zp = double(f(xp,yp));
for i = 2:kterm
    xp_n = xp;
    yp_n = yp;
    zp_n = zp;
    v = lista_punti{i,1};
    xp = v(1);
    yp = v(2);
    zp = double(f(xp,yp));
    pl = line([xp_n xp],[yp_n yp],[zp_n zp]);   %visualizzo le iterazioni dell'algoritmo
    pl.Marker = '*';
    pl.Color = 'red';
    pl.LineWidth = 1;
end
surf(X,Y,Z);

xp = xk(1);
yp = xk(2);
zp = double(f(xp,yp));
plot3(xp,yp,zp,'y*');
legend(strcat('cond  = ',int2str(mu)),strcat('iter = ',int2str(kterm)));
view(3);
hold off


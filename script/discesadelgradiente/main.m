import discesa.*


N = 100;
mu = 100000; %indice di condizionamento
A = full(sprandsym(N, 1, 1/mu, 1)) * 100;  %costruisco la matrice simmetrica e definita positiva 
                                      %(dim, densità, 1/indice_condizionamento, definita positiva = 1) 
               
%parametri
b = rand(N,1) * 100;
x0 = rand(N,1) * 100;      %oppure considero il vettore nullo come pos iniziale
nmax = 1000;
toll = 10^(-8);

%variabili per memorizzare dati
kterm = 0;
lista_punti = cell(nmax,1);


%algoritmo
[xk,lista_punti,kterm] = discesa(A, b, x0, nmax, toll,lista_punti);

%valutazione errore rispetto al metodo diretto
xt = A\b;
ea = norm(xk-xt);
er = norm(xk-xt)/norm(xt);
rk = b-A*xk;


%Andamento del residuo relativo norm(r)/norm(b)
res_vect = zeros(kterm,1);
for i = 1:kterm
    xc = lista_punti{i,1};
    rc = (norm(b-A*xc)/norm(b));
    res_vect(i) = rc;
end


%grafico del residuo
figure(1);
hold on
semilogy(1:kterm,res_vect,'blue');
xlabel('Numero di iterazioni')
ylabel('Residuo relativo')
legend(strcat('cond = ',int2str(mu)))
set(gca, 'YScale', 'log')
hold off


import gradienteconiugato.*
import p_gradienteconiugato.*

N = 100;
mu = 400; %indice di condizionamento 1
A = full(sprandsym(N, 1, 1/mu, 1)) * 100;  %costruisco la matrice simmetrica e definita positiva 
                                      %(dim, densità, 1/indice_condizionamento, definita positiva = 1) 
               
%parametri
b = rand(N,1) * 100;
x0 = rand(N,1) * 100;      %oppure considero il vettore nullo come pos iniziale
nmax = 1000;
toll = eps(norm(b));

%variabili per memorizzare dati
kterm = 0;
lista_punti = cell(nmax,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%algoritmo
[xk,lista_punti,kterm] = gradienteconiugato(A, b, x0, nmax, toll,lista_punti);

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = max(sum(abs(A),2)./diag(A))-2; %per evitare errore: Encountered nonpositive pivot. all'alumentare di mu
R1 = ichol(sparse(A), struct('type','ict','droptol',1e-3,'diagcomp',alpha)); %fatt cholesky incompleta
R2 = R1';
%algoritmo
[xk2,lista_punti2,kterm2] = p_gradienteconiugato(A, b, x0, nmax, toll,R1,R2,lista_punti);

%valutazione errore rispetto al metodo diretto
ea2 = norm(xk2-xt);
er2 = norm(xk2-xt)/norm(xt);
rk2 = b-A*xk2;


%Andamento del residuo relativo norm(r)/norm(b)
res_vect2 = zeros(kterm2,1);
for i = 1:kterm2
xc2 = lista_punti2{i,1};
rc2 = (norm(b-A*xc2)/norm(b));
res_vect2(i) = rc2;
end


%grafico del residuo
legend ;
hold on
semilogy(1:kterm,res_vect,'blue');

semilogy(1:kterm2,res_vect2,'red');

legend('no precond','precond')



xlabel('Numero di iterazioni')
ylabel('Residuo relativo')
plot(nan, nan, 'DisplayName', strcat('cond = ',int2str(mu)))
set(gca, 'YScale', 'log')
hold off


import gradienteconiugato.*
import p_gradienteconiugato.*
import discesa.*

N = 100;
mu = 5000; %indice di condizionamento 1
A = full(sprandsym(N, 1, 1/mu, 1)) * 100;  %costruisco la matrice simmetrica e definita positiva 
                                      %(dim, densità, 1/indice_condizionamento, definita positiva = 1) 

%parametri
b = rand(N,1) * 100;
x0 = rand(N,1) * 100;      %oppure considero il vettore nullo come pos iniziale
nmax = 1000;
toll = 1e-12;
toll2 = toll * norm(b);

%variabili per memorizzare dati
lista_punti = cell(nmax,1);

%metodo diretto
xt = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%algoritmo CG
[xk,lista_punti,kterm] = gradienteconiugato(A, b, x0, nmax, toll2,lista_punti);

    %Andamento del residuo relativo norm(r)/norm(b)
res_vect = zeros(kterm,1);
for i = 1:kterm
xc = lista_punti{i,1};
rc = (norm(b-A*xc)/norm(b));
res_vect(i) = rc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%algoritmo PCG
alpha = max(sum(abs(A),2)./diag(A))-2; %per evitare errore: Encountered nonpositive pivot.
R1 = ichol(sparse(A), struct('type','ict','droptol',1e-3,'diagcomp',alpha)); %fatt cholesky incompleta
R2 = R1';

[xk2,lista_punti2,kterm2] = p_gradienteconiugato(A, b, x0, nmax, toll2,R1,R2,lista_punti);

    %Andamento del residuo relativo norm(r)/norm(b)
res_vect2 = zeros(kterm2,1);
for i = 1:kterm2
xc = lista_punti2{i,1};
rc = (norm(b-A*xc)/norm(b));
res_vect2(i) = rc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%algoritmo Discesa del gradiente
[xk3,lista_punti3,kterm3] = discesa(A, b, x0, nmax, toll2,lista_punti);


%Andamento del residuo relativo norm(r)/norm(b)
res_vect3 = zeros(kterm3,1);
for i = 1:kterm3
xc = lista_punti3{i,1};
rc = (norm(b-A*xc)/norm(b));
res_vect3(i) = rc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%algoritmo PCG matlab

[xk4,flag,rel_res,kterm4,res_vect4] = pcg(A, b, toll, nmax,R1,R2,x0);
kterm4 = kterm4+1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%algoritmo gmres matlab

[xk5,flag,rel_res,kterm5,res_vect5] = gmres(A, b,N, toll, size(A,1),R1,R2,x0);
kterm5 = kterm5+1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
semilogy(1:kterm,res_vect,'blue','LineWidth', 1.5);
semilogy(1:kterm2,res_vect2,'red','LineWidth', 1.5);
semilogy(1:kterm3,res_vect3,'black','LineWidth', 1.5);
semilogy(1:kterm4,res_vect4,'green','LineWidth', 1.5);
semilogy(1:kterm5(2),res_vect5,'cyan','LineWidth', 1.5);

legend('cg','pcg','discesa','pcg matlab','gmres matlab')


xlabel('Numero di iterazioni')
ylabel('Residuo relativo')
plot(nan, nan, 'DisplayName', strcat('cond = ',int2str(mu)))
set(gca, 'YScale', 'log')

hold off

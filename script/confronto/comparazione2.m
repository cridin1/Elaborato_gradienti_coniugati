import gradienteconiugato.*
import p_gradienteconiugato.*
import discesa.*

format short
mu = [5 50 500 5000 50000]'; %indice di condizionamento 1
ea = zeros(5,5);
er = zeros(5,5);
kmat = zeros(5,5);
tempo = zeros(5,5);

for i = 1:size(mu,1)
    N = 100;
    A = full(sprandsym(N, 1, 1/mu(i), 1)) * 100;  %costruisco la matrice simmetrica e definita positiva 
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
    tic;
    [xk,lista_punti,kterm] = gradienteconiugato(A, b, x0, nmax, toll2,lista_punti);
    tempo(i,1) = toc;
    kmat(i,1) = kterm;
    ea(i,1) = norm(xk-xt);
    er(i,1) = ea(i,1)/norm(xt);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %algoritmo PCG
    tic;
    alpha = max(sum(abs(A),2)./diag(A))-2; %per evitare errore: Encountered nonpositive pivot.
    R1 = ichol(sparse(A), struct('type','ict','droptol',1e-3,'diagcomp',alpha)); %fatt cholesky incompleta
    R2 = R1';

    [xk2,lista_punti2,kterm2] = p_gradienteconiugato(A, b, x0, nmax, toll2,R1,R2,lista_punti);

    tempo(i,2) = toc;
    kmat(i,2) = kterm2;
    ea(i,2) = norm(xk2-xt);
    er(i,2) = ea(i,2)/norm(xt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %algoritmo Discesa del gradiente
    tic;
    [xk3,lista_punti3,kterm3] = discesa(A, b, x0, nmax, toll2,lista_punti);
    tempo(i,3) = toc;
    kmat(i,3) = kterm3;
    ea(i,3) = norm(xk3-xt);
    er(i,3) = ea(i,3)/norm(xt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %algoritmo PCG matlab
    tic;
    [xk4,flag,rel_res,kterm4,res_vect4] = pcg(A, b, toll, nmax,R1,R2,x0);
    kterm4 = kterm4+1;
    tempo(i,4) = toc;
    kmat(i,4) = kterm4;
    ea(i,4) = norm(xk4-xt);
    er(i,4) = ea(i,4)/norm(xt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %algoritmo gmres matlab
    tic;
    [xk5,flag,rel_res,kterm5,res_vect5] = gmres(A, b,N, toll, size(A,1),R1,R2,x0);
    kterm5 = kterm5+1;
    tempo(i,5) = toc;
    kmat(i,5) = kterm5(2);
    ea(i,5) = norm(xk5-xt);
    er(i,5) = ea(i,5)/norm(xt);
    
end

ea = round(ea,4,"significant");
er =round(er,4,"significant");
kmat = round(kmat,4,"significant");
tempo = round(tempo,4,"significant");

ea = [mu ea];
er = [mu er];
kmat = [mu kmat];
tempo = [mu tempo];

writematrix(ea,'ea.csv')
writematrix(er,'er.csv')
writematrix(kmat,'kmat.csv')
writematrix(tempo,'tempo.csv')


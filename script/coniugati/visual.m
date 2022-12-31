import gradienteconiugato.*
import p_gradienteconiugato.*
import discesa.*


N = 2;
mu = 5;
A = full(sprandsym(N, 0.8, 1/mu, 1))*10;  %costruisco la matrice simmetrica e definita positiva 
                                      %(dim, densità, 1/indice_condizionamento, definita positiva = 1) 


%parametri
b = rand(N,1)*10;
x0 = rand(N,1)*10;       %considero il vettore nullo come pos iniziale
nmax = 80;
toll = 1e-12;

lista_punti = cell(nmax,1);

%met. diretto
xt = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Alg1 cg
%algoritmo
[xk,lista_punti,kterm] = gradienteconiugato(A, b, x0, nmax, toll,lista_punti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Alg2 pcg

alpha = 0;%max(sum(abs(A),2)./diag(A))-2; %per evitare errore: Encountered nonpositive pivot.
R1 = ichol(sparse(A));%, struct('type','ict','droptol',1e-3,'diagcomp',alpha)); %fatt cholesky incompleta
R2 = R1';
%algoritmo
[xk2,lista_punti2,kterm2] = p_gradienteconiugato(A, b, x0, nmax, toll,R1,R2,lista_punti);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Alg3 discesa
%algoritmo
[xk3,lista_punti3,kterm3] = discesa(A, b, x0, nmax, toll,lista_punti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%visualizzazione
syms x y
f = symfun(0.5*(A(1,1) * x^2 + (A(1,2)+A(2,1)) * x*y + A(2,2) * y^2) - b(1)*x - b(2)*y,[x y]);

%tuning dei parametri per la visualizzazione
camp = 0.5; %passo di campionamento
off = norm(x0); %scelgo offset
[X,Y] = meshgrid(-off+xt(1):camp:off+xt(1),-off+xt(2):camp:off+xt(2)); %griglia di punti


Z = 0.5*(A(1,1) * X.*X + (A(1,2)+A(2,1)) * X.*Y + A(2,2) * Y.*Y) - b(1)*X - b(2)*Y;
figure(1);
hold on

%%%%Cg
v = lista_punti{1,1};
xp = v(1);
yp = v(2);
zp = double(f(xp,yp));
for i = 2:(kterm)
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
    pl.LineWidth = 2;
end

% %pcg
% v = lista_punti2{1,1};
% xp = v(1);
% yp = v(2);
% zp = double(f(xp,yp));
% for i = 2:(kterm2)
%     xp_n = xp;
%     yp_n = yp;
%     zp_n = zp;
%     v = lista_punti2{i,1};
%     xp = v(1);
%     yp = v(2);
%     zp = double(f(xp,yp));
%     pl2 = line([xp_n xp],[yp_n yp],[zp_n zp]);   %visualizzo le iterazioni dell'algoritmo
%     pl2.Marker = '*';
%     pl2.Color = 'green';
%     pl2.LineWidth = 2;
% end

%%Discesa
v = lista_punti3{1,1};
xp = v(1);
yp = v(2);
zp = double(f(xp,yp));
for i = 2:kterm3
    xp_n = xp;
    yp_n = yp;
    zp_n = zp;
    v = lista_punti3{i,1};
    xp = v(1);
    yp = v(2);
    zp = double(f(xp,yp));
    pl3 = line([xp_n xp],[yp_n yp],[zp_n zp]);   %visualizzo le iterazioni dell'algoritmo
    pl3.Marker = '*';
    pl3.Color = 'cyan';
    pl3.LineWidth = 1;
end


surf(X,Y,Z);

xp = xt(1);
yp = xt(2);
zp = double(f(xp,yp));
plot3(xp,yp,zp,'y*');
plot(nan,nan,'y');

view(3);
hold off
legend(strcat('iter cg = ',int2str(kterm)),strcat('iter pcg = ',int2str(kterm2)),strcat('iter discesa = ',int2str(kterm3)),strcat('cond = ',int2str(mu)));


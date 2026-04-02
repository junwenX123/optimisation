clear
close all


methode="SR1";

%%% Test 1 : fonctionnelle quadratique (1/2)<Qx,Qx> + <g,x> sous la 
%%% contrainte Ax-b=0. 

Q=eye(2,2);
g=[0;0];
A=[1,1];
b=1;
f1=@(x) 0.5*(Q*x).'*(Q*x) + g.'*x;
c1=@(x) A*x-b;
F1=@(x) [f1(x); c1(x)];

x=[5;4];
tol=10^(-6);
pas_gradient=[0.001;0.001];
max_iter=100;
nfonc_max=250;

[~, ~, tableau1]=SQP(F1,x,tol, pas_gradient ,max_iter,nfonc_max,methode);
display(tableau1);


%%% Test 2 : MHW4D

f2=@(x) (x(1)-1)^2 + (x(1)-x(2))^2+(x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;
c2=@(x) [ x(1)+x(2)^2+x(3)^2-3*sqrt(2)-2 ; x(2)-x(3)^2+x(4)-2*sqrt(2)+2 ; x(1)*x(5)-2 ];
F2=@(x) [f2(x); c2(x)];

x=[-1;2;1;-2;-2];
tol=10^(-6);
pas_gradient=[0.001;0.001;0.001;0.001;0.001];
max_iter=100;
nfonc_max=250;
borne_inf=[-2;1;0;-3;-3];
borne_sup=[0;3;2;0;-1];

[~, ~, tableau2]=SQP(F2,x,tol, pas_gradient, max_iter,nfonc_max,methode,borne_inf,borne_sup);
display(tableau2);


%%% Test 3 : Ariane1

mu=1700;
k=[0.1101; 0.1532; 0.2154];
ve=[2647.2; 2922.4; 4344.3];
delta_requis=11527;
F3=@(x) f_ariane1(x, mu, k, ve, delta_requis);

x_sol=[145349;31215; 7933];
x1=[200000;50000; 10000];
x2=[140000;30000; 7000];
x3=[10000; 5000; 1000];
tol=10^(-6);
pas_gradient=[1;1;1];
max_iter=100;
nfonc_max=350;
borne_inf=[0;0;0];
borne_sup=[inf;inf;inf];

[~, ~, tableau3]=SQP(F3,x1,tol, pas_gradient, max_iter, nfonc_max, methode,borne_inf,borne_sup);
display(tableau3);

%table2latex(tableau1, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\BFGS_quadra');
%table2latex(tableau2, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\BFGS_MWH4D');
%table2latex(tableau3, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\BFGS_Ariane1');
%fclose('all');
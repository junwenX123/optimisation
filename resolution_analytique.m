clear
close all

currentFile=mfilename("fullpath");
globalFolder=fileparts(currentFile);
simulateurFolder=strcat(globalFolder,"/Simulateur_trajectoire");
SQPFolder=strcat(globalFolder,"/SQP");
addpath(globalFolder);
addpath(simulateurFolder);
addpath(SQPFolder);

%%% Initialisation des données 

k=[0.10; 0.15; 0.20];
ve=[2600; 3000; 4300];
mu=3.986*10^(14);
Rt=6378137;
Rc=Rt+250000;
vc=sqrt(mu/Rc);
vp=1.2*vc;

%%% Définition de la fonction sur laquelle les itérations de Newton seront
%%% appliquées 

F_newton=@(x3) f_newton(x3, k, ve, vp);

%%% Itérations de Newton 

x3_0=2.5;
tol=10^(-8);
max_iter=100;
pas_gradient=10^(-4);

[x3, c_x3, iter] = newton_ralphson(F_newton, x3_0, tol, max_iter, pas_gradient);

%%% Calculs de x1, x2, les yi et lambda

Omega = k ./ (1 + k);  
x1 = (1 - (ve(3)/ve(1)) * (1 - Omega(3)*x3)) / Omega(1);
x2 = (1 - (ve(3)/ve(2)) * (1 - Omega(3)*x3)) / Omega(2);
x=[x1;x2;x3];
y=(1+k)./x-k;
f_x=-y(1)*y(2)*y(3);
lambda = f_x/(ve(3)*(1 - Omega(3)*x(3)));


%%% Calcul des composantes du gradient du lagrangien (trois premières
%%% égalités des conditions KKT) et vérification de l'égalité à 0

grad_l=[1;1;1];

for j=1:3
    product_y=1;
    for i=1:3
        if i~=j
            product_y=product_y*y(i);
        end
    end
    grad_l(j)= product_y*(1+k(j))/(x(j)^2);
    grad_l(j)=grad_l(j)+ lambda*ve(j)/x(j);
end

%%% Affichage des données obtenues 

display(x1);
display(x2);
display(x3);
display(lambda);
display(c_x3);
display(iter);
display(grad_l);

%%% Récupération des masses d'ergols

masse_u=1000;
[me,Mtot]=me_from_x(x, k, masse_u);
display(me);
display(Mtot);
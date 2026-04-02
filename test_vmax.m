clear
close all

currentFile=mfilename("fullpath");
globalFolder=fileparts(currentFile);
simulateurFolder=strcat(globalFolder,"/Simulateur_trajectoire");
SQPFolder=strcat(globalFolder,"/SQP");
addpath(globalFolder);
addpath(simulateurFolder);
addpath(SQPFolder);

%%% Déterminer les masses d'ergol pour notre problème 
%%% masse_utile=1000kg
%%% altitude de l'orbite=250kg

methode="BFGS";
masse_u=1000;
k=[0.10; 0.15; 0.20];
ve=[2600; 3000; 4300];
mu=3.986*10^(14);
Rt=6378137;
Rc=Rt+250000;
vc=sqrt(mu/Rc);
vp=1.2*vc;
F_etagement=@(x) -f_etagement(x, masse_u, k, ve, vp);

x=[19000; 7000; 2000];
tol=10^(-8);
pas_gradient=[0.001;0.001;0.001];
max_iter=100;
nfonc_max=250;
borne_inf=[0;0;0];
borne_sup=[inf;inf;inf];

[sol, F_sol, tableau_etagement]=SQP(F_etagement,x, tol, pas_gradient, max_iter, nfonc_max, methode,borne_inf,borne_sup);
display(tableau_etagement);

me=sol;
masse_lanceur=-(F_sol(1)/masse_u)^(-1);

%%% Estimation des paramètres angulaires pour la première résolution de
%%% trajetoire

% Angles convenables déterminés manuellement

theta0=30*pi/180;
theta1=-10*pi/180;
theta2=-37*pi/180;
theta3=-26*pi/180;

theta=[theta0; theta1; theta2; theta3];

% Vitesse, trajectoire et masse initiales

V0=100.*[cos(theta0);sin(theta0)];
R0= [6378137;0];
M0=masse_lanceur;

% Autres donnees

ve=[2600; 3000; 4400]; % vitesses d'éjection
alpha=[15; 10; 10]; % accélérations à l'allumage
k=[0.10; 0.15; 0.20]; % indices constructifs

% Mise en forme pour le simulateur

donnee=zeros(5,3);
donnee(1,:)=me;
donnee(2,:)=ve;
donnee(3,:)=alpha;
donnee(4,:)=theta(2:end);
donnee(5,:)=k;

% Calcul de cette première trajectoire

[time, result, tf, Mf, Rf, Vf]=simulateur_trajectoire(R0, V0, M0, donnee);
tracer_traj(time, result, Mf);

%%% Optimisation de la vitesse

facteur_echelle=10^(-3);
option=odeset(RelTol=1e-6,AbsTol=1e-8);

F_vitesse=@(x) -facteur_echelle*f_vitesse(x, R0, M0, donnee, Rc, option);

theta=[theta0; theta1; theta2; theta3];
tol=10^(-6);
pas_gradient=[10^(-6);10^(-6);10^(-6);10^(-6)];
max_iter=100;
nfonc_max=300;
methode="BFGS";


[theta_opti, v_opti_moins, tableau_vitesse]=SQP(F_vitesse,theta, tol, pas_gradient, max_iter, nfonc_max, methode);
v_opti=-v_opti_moins(1);
display(tableau_vitesse);

%table2latex(tableau_etagement, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\etagement_iter1_vmax')
%table2latex(tableau_vitesse, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\vitesse_iter1_v_max')


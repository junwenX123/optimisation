clear
close all

currentFile=mfilename("fullpath");
globalFolder=fileparts(currentFile);
simulateurFolder=strcat(globalFolder,"/Simulateur_trajectoire");
SQPFolder=strcat(globalFolder,"/SQP");
addpath(globalFolder);
addpath(simulateurFolder);
addpath(SQPFolder);

%%% Paramètres globaux

% Cahier des charges

masse_u=1000;
R_orbite=250000;

% Etagement

methode_etagement="BFGS";
k=[0.10; 0.15; 0.20];
ve=[2600; 3000; 4300];
mu=3.986*10^(14);
Rt=6378137;
Rc=Rt+R_orbite;
vc=sqrt(mu/Rc);
vp=1.2*vc;

x=[19000; 7000; 2000];
tol=10^(-6);
pas_gradient_etagement=[0.001;0.001;0.001];
max_iter=100;
nfonc_max=350;
borne_inf=[0;0;0];
borne_sup=[inf;inf;inf];

% Maximisation de la vitesse

methode_vitesse="BFGS";
R0=[6378137;0];
alpha=[15; 10; 10];
pas_gradient_vitesse=[10^(-6);10^(-6);10^(-6);10^(-6)];
option=odeset(RelTol=1e-6,AbsTol=1e-8);

% Angles de départ expérimentaux pour le lanceur à vp=1.2vc

theta0=30*pi/180;
theta1=-10*pi/180;
theta2=-37*pi/180;
theta3=-26*pi/180;
theta=[theta0; theta1; theta2; theta3];

% Mise en forme des données pour les optimisations de trajectoires

donnee=zeros(5,3);
donnee(2,:)=ve;
donnee(3,:)=alpha;
donnee(4,:)=theta(2:end);
donnee(5,:)=k;


%%% Recherche du lanceur optimal ave un facteur d'échelle

iter=1;
delta_v=100;
facteur_echelle=10^(-3);

while abs(delta_v) > 10^(-4) && iter<100
    
    if iter==1
        delta_v=0;
    end
    vp=vp+delta_v;

    fprintf("Résolution du problàme !\n");
    fprintf("Iteration : %d\n", iter);
    fprintf("Facteur d'échelle :  %f\n", facteur_echelle)
    fprintf("Vitesse propulsive :  %f\n", vp)


    % Etagement
    
    F_etagement=@(x) -f_etagement(x, masse_u, k, ve, vp);
    [sol, F_sol, tableau_etagement]=SQP(F_etagement,x, tol, pas_gradient_etagement, max_iter, nfonc_max, methode_etagement,borne_inf,borne_sup);
    me=sol;
    masse_lanceur=-(F_sol(1)/masse_u)^(-1);
    donnee(1,:)=me;
    %display(masse_lanceur)
    display(tableau_etagement) 
    %table2latex(tableau_etagement, strcat('C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\etagement_iter',sprintf('%d',iter)))

    % Vitesse maximale
   
    F_vitesse_echelle=@(x) -facteur_echelle.*f_vitesse(x, R0, masse_lanceur, donnee, Rc,option);
    [theta_opti, v_opti_moins, tableau_vitesse]=SQP(F_vitesse_echelle,theta, tol, pas_gradient_vitesse, max_iter, nfonc_max, methode_vitesse);
    v_opti=-v_opti_moins(1)/facteur_echelle;
    display(tableau_vitesse)
    %table2latex(tableau_vitesse, strcat('C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\vitesse_iter',sprintf('%d',iter)))
    
    % Changement de la vitesse de propulsion et des angles initiaux
    
    theta=theta_opti;
    delta_v=vc-v_opti;
    iter=iter+1;
end


% Recherche des angles optimaux sans le facteur d'échelle

fprintf("Résolution finale sans facteur d'échelle !\n");
fprintf("Vitesse propulsive :  %f\n", vp)

donnee(4,:)=theta(2:end);
option=odeset(RelTol=1e-6,AbsTol=1e-8);
F_vitesse=@(x) -f_vitesse(x, R0, masse_lanceur, donnee, Rc, option);
pas_gradient_vitesse=[10^(-5);10^(-5);10^(-5);10^(-5)];
[theta_opti, v_opti_moins, tableau_vitesse]=SQP(F_vitesse,theta, tol, pas_gradient_vitesse, max_iter, nfonc_max, methode_vitesse);
v_opti=-v_opti_moins(1);
display(tableau_vitesse)
%table2latex(tableau_vitesse, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\vitesse_finale')

% Simulation de la trajectoire avec les angles optimaux sans facteur
% d'échelle

donnee(4,:)=theta_opti(2:end);
V0=100.*[cos(theta_opti(1));sin(theta_opti(1))];
[time, result, tf, Mf, Rf, Vf]=simulateur_trajectoire(R0, V0, masse_lanceur, donnee, option);
altitude_finale= norm(Rf, 2);
vitesse_finale= norm(Vf, 2);

%%% Tableau des caractéristiques du lanceur optimisé et graphiques de la
%%% trajectoire

% Réecriture des données vectorielles dans une chaîne de caractères afin
% de pouvoir les afficher dans un tableau

me_str="(";

for i=1:length(me)
    if i<length(me)
        me_str=strcat(me_str,sprintf("%f;",me(i)));
    else
        me_str=strcat(me_str,sprintf("%f)",me(i)));
    end
end

theta_opti=theta_opti.*(180/pi);
theta_opti_str="(";

for i=1:length(theta_opti)
    if i<length(theta_opti)
        theta_opti_str=strcat(theta_opti_str,sprintf("%f;",theta_opti(i)));
    else
        theta_opti_str=strcat(theta_opti_str,sprintf("%f)",theta_opti(i)));
    end
end

% Fabrication des tableaux

cahier_des_charges=table(masse_u, R_orbite*10^(-3), vc, 'VariableNames', ...
                    ["Masse utile (en kg)", "Altitude de l'orbite (en km)", "Vitesse cible (en m/s)"]);
caracteristiques_etagement=table(vp, me_str, masse_lanceur, 100*masse_u/masse_lanceur, 'VariableNames', ...
                     ["Vitesse propulsive (en m/s)", "Masses d'ergol par étage (en kg)", "Masse totale du lanceur (en kg)","Pourcentage de masse utile"]);

var_names=[ "Angles optimaux (en degré)","Masse finale (en kg) ", ...
        "Erreur sur l'altitude (en m)", "Erreur sur la vitesse (en m/s)", ...
        "Pr. sca. position et vitesse finales (en m^2/s)"];
caracteristiques_trajectoire=table( theta_opti_str, Mf, altitude_finale-Rc, vitesse_finale-vc, ...
                 Rf.'*Vf,'VariableNames', var_names);

% Affichages 

display(cahier_des_charges);
display(caracteristiques_etagement);
display(caracteristiques_trajectoire);
tracer_traj(time, result, Mf);

%table2latex(cahier_des_charges, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\cahier_charge')
%table2latex(caracteristiques_etagement, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\caracteristiques_etagement')
%table2latex(caracteristiques_trajectoire, 'C:\Users\minib\Documents\Cours\M2_reprise_etudes\Bloc_2\Optimisation\Projet optimisation\Rapport\caracteristiques_trajectoire')
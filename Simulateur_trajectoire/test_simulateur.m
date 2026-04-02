clear
close all


%%% Données initiales

% Vitesse, trajectoire et masse

V0=100.*[cos(15);sin(15)];
R0= [6378137;0];
M0=208611;

% Ariane1

%%% Mise en forme des données pour le simulateur

donnee=zeros(5,3);

% Masses d'ergol par étage

donnee(1,1)=145349;
donnee(1,2)=31215;
donnee(1,3)=7933;

% Vitesses d'éjection par étage

donnee(2,1)=2647.2;
donnee(2,2)=2922.4;
donnee(2,3)=4344.3;

% Accélération à l'allumage par étage

donnee(3,1)=15;
donnee(3,2)=10;
donnee(3,3)=10;

% Angles de poussée par étage

donnee(4,1)=25*pi/180;
donnee(4,2)=-60*pi/180;
donnee(4,3)=-31*pi/180;

% Indice constructif par étage

donnee(5,1)=0.1101;
donnee(5,2)=0.1532;
donnee(5,3)=0.2154;

%%% Calculs

[time, result, tf, Mf, Rf, Vf]=simulateur_trajectoire(R0, V0, M0, donnee);

%%% Affichage des graphiques

tracer_traj(time, result, Mf);

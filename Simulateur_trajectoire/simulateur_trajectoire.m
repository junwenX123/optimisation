function [time, result, tf, Mf, Rf, Vf]=simulateur_trajectoire(R0, V0, M0, donnee, option)
% R0 la position initiale
% V0 la vitesse initiale
% M0 la masse initiale
% donnee le tableau (5,3) contenant à la colonne j les données de l'étage j
% dans l'ordre suivant : masse d'ergol de l'étage j, vitesse d'éjection de
% l'étage j, l'accélération demandée à l'allumage de l'étage j, l'angle
% d'incidence thetaj, kj l'indice constructif de l'étage j. 
% option les paramètres d'options à utiliser pour ode45
% time le vecteur des instants utilisés pour résoudre l'EDO.
% result la matrice dont les lignes contiennent à l'instant correspondant
% les valeurs de [R, V, M] (la dernière valeurde M est la masse avant
% décrochement du troisième étage).
% tf l'instant final, Mf la masse finale une fois le troisième étage
% détaché, Rf l'altitude finale et Vf la vitesse finale

Ri=R0;
Vi=V0;
Mi=M0;
ti=0;
time=[];
result=[];

for j=1:3


    % Récupération des données de l'étage j

    mej=donnee(1,j);
    vj=donnee(2,j);
    aj=donnee(3,j);
    thetaj=donnee(4,j);
    kj=donnee(5,j);
    Tj=aj*Mi;
    qj=Tj/vj;
 

    % Création de la fonction à donner en argument à ode_traj pour 
    % l'étage j
    
    ode_fun=@(t, y) f_traj(t, y, Tj, qj, thetaj);

    % Résolution de l'EDO

    y0=[Ri; Vi; Mi];
    tf=ti+mej/qj;
    if nargin==4
        [t_span, yf]=ode45(ode_fun,[ti; tf], y0);
    else 
        [t_span, yf]=ode45(ode_fun,[ti; tf], y0, option);
    end
    

    % Récupération des données 

    time=[time; t_span];
    result=[result; yf];


    % Séparation de la masse structurelle de l'étage et initialisation
    % des variables

    ti=t_span(end);
    Ri=yf(end,1:2).';
    Vi=yf(end,3:4).';
    Mi=yf(end,5);
    Mi=Mi-kj*mej;

    
end

tf=ti;
Rf=Ri;
Vf=Vi;
Mf=Mi;
end

function [c]=f_newton(x3, k, ve, vp)
% x3 la variable de la résolution analytique telle qu'obtenue dans le
% rapport du projet.
% k le vecteur colonne des indices constructifs du lanceur.
% ve le vecteur colonne des vitesses d'éjection du lanceur.
% vp la vitesse propulsive du lanceur.
% c les contraintes du système KKT reformulées en fonction de x3. 

Omega = k ./ (1 + k);  
x1 = (1 - (ve(3)/ve(1)) * (1 - Omega(3)*x3)) / Omega(1);
x2 = (1 - (ve(3)/ve(2)) * (1 - Omega(3)*x3)) / Omega(2);
x=[x1;x2;x3];

% Calcul des contraintes c(x1,x2,x3)

c=0;
for j=1:3
    c=c+ve(j)*log(x(j));
end
c=c-vp;
end
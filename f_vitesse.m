function [F_theta]=f_vitesse(theta, R0, M0, donnee, Rc, option)
% theta les angles de poussée
% R0, M0, donnee les arguments initiaux de simulateur_trajectoire
% Rc l'altitude sible
% option les paramètres d'options à utiliser pour ode45
% F_theta le vecteur colonne correspondant à [fonction; contrainte] pour 
% le problème de maximisation de la vitesse. 

donnee(4,:)=theta(2:end);
V0=100.*[cos(theta(1));sin(theta(1))];
if nargin==5
    [~, ~, ~, ~, Rtf, Vtf]=simulateur_trajectoire(R0, V0, M0, donnee);
else 
    [~, ~, ~, ~, Rtf, Vtf]=simulateur_trajectoire(R0, V0, M0, donnee, option);
end
F_theta=[norm(Vtf,2); norm(Rtf,2)-Rc; (Rtf.')*Vtf];
end
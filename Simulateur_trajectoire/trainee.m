function [D]=trainee(rho, V)
% rho la densité de l'atmosphère
% V le vecteur vitesse
% D le vecteur de traînée

cx=0.1; % coefficient de trainée du lanceur
D=(-cx*rho*norm(V,2)).*V;
end
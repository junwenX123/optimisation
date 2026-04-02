function [ode_vec]=f_traj(t, y, Tj, qj, thetaj)
% t l'instant considéré, y le vecteur [R; V; M] contenant l'altitude,
% la vitesse et la masse courantes.
% Tj la poussée de l'étage j, vej sa vitesse d'éjection, qj le débit
% associé et thetaj son angle d'incidence.
% ode_vec le vecteur correspondant au membre de droite de l'ODE.


% Récupération de l'altitude, vitesse et masse courantes 

R=y(1:2);
V=y(3:4);
M=y(5);

% Calcul du membre de doite de l'ODE du problème.

rho=atmosphere_density(R);
W=weight(R, M);
D=trainee(rho, V);
T=poussee(R, V, Tj, thetaj);
ode_vec=zeros(5,1);
ode_vec(1:2,1)=V; 
ode_vec(3:4,1)=(T+W+D)./M;
ode_vec(5,1)=-qj;

end

    
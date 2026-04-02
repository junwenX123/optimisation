function [rho]=atmosphere_density(R)
% R la position
% rho la dentisté de l'atmosphère à l'altitude norm(R,2)

rho0=1.225; % densité au sol en kg/m³
Rt=6378137; % rayon terrestre en m
H=7000; % facteur d'écchelle
altitude=norm(R, 2);
rho=rho0*exp(-(altitude-Rt)/H);
end
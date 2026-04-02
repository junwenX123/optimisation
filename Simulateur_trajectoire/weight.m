function [W]=weight(R, M)
% R la position
% M la masse 
% W le vecteur de poids

mu=3.9861014; % constante gravitationnelle terrestre en m³/s²
altitude=norm(R,2);
W=(-mu*M/altitude^3).*R;
end
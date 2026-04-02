function [T]=poussee(R, V, Tj, thetaj)
% R la position
% V la vitesse
% Tj la poussée de l'étage j
% thetaj l'angle d'incidence
% T le vecteur de poussée

altitude=norm(R,2);
gamma=asin(R.'*V/(altitude*norm(V,2)));
er=(1/altitude).*R;
eh=(1/altitude).*[-R(2);R(1)];
u=cos(gamma+thetaj).*eh + sin(gamma+thetaj).*er;
T=Tj.*u;
end
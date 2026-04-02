function [new_H]=modif_defpos(H,epsilon)
% H la matrice symétrique qu'on veut rendre définie positive. 
% epsilon paramètre de modification strictement positifs des valeurs propres.
% new_H la matrice symétrique modifiée proche de H et définie positive.
if nargin==1
    epsilon=max(H)*10^(-1);
end
vp=eig(H);
min_vp=min(vp);
if min_vp <= 0
    new_H=H+(abs(min_vp)+epsilon)*eye(size(H,1));
else
    new_H=H;
end
end
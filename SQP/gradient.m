function [grad]=gradient(f, x, f_x, delta)
% f la fonction dont on calcule le gradient.
% x le point où l'on calcule le gradient (vecteur colonne).
% delta le vecteur colonne des déplacements pour l'approximation de type 
% différences finies pour chaque variable.
% f_x le vecteur (colonne) f(x).
% delta et f_x sont des arguments optionnels.
% grad la matrice des gradients (colonnes) au point x (donc la
% transposée de la jacobienne).

pas=zeros(length(x),1);
grad=zeros(length(x),length(f_x));
for i=1:length(x)
    pas(i)=delta(i);
    grad(i,:)=(f(x+pas)-f_x)./delta(i);
    pas(i,:)=0;
end
end
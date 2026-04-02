function [x_pas, nfonc]=armijo(F, merite, rho, d, x, merite_x, grad_merite_x,nfonc, borne_inf,borne_sup)
% F la fonction de SQP
% merite la fonction mérite.
% rho le paramètre de pénalisation.
% d une direction de descente au point x. 
% merite_x la fonction mérite calculée au point x.
% grad_merite_x la dérivée partielle de la fonction mérite selon d au point
% x. 
% borne_inf les bornes inférieures sur la variable, borne_sup ses bornes
% supérieures. Les deux arguments sont optionnels si il n'y a pas 
% de borne sur la variable mais les deux doivent l'être sinon.
% x_pas le prochain point atteint avec en ajoutant pas.


% Initialisation

pas=1;
coef=0.1;
x_pas=x+pas*d;
F_pas=F(x_pas);
nfonc=nfonc+1;
norme_pas=norm(F_pas(2:end), 1);
merite_pas=merite(rho,F_pas,norme_pas); 
iter=0;


% Itérations d'Armijo

while merite_pas > merite_x+coef*pas*grad_merite_x && iter<=10
    iter=iter+1;
    pas=pas/2;
    x_pas=x+pas*d;
    F_pas=F(x_pas);
    nfonc=nfonc+1;
    norme_pas=norm(F_pas(2:end), 1);
    merite_pas=merite(rho,F_pas,norme_pas); 
end

% Projection si il y a des bornes

if nargin > 8
    for i=1:length(x_pas)
        x_pas(i)=max(x_pas(i),borne_inf(i));
        x_pas(i)=min(x_pas(i),borne_sup(i));
    end
end

end
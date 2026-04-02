function [sol, F_sol, tableau]=SQP(F, x, tol, pas_gradient, max_iter, max_nfonc, methode, borne_inf, borne_sup)
% F=[f;c] où f est la fonction à minimiser, c sont les contraintes d'égalité. 
% x le point de départ.
% tol la tolerance sur la norme du gradient du lagrangien. 
% pas_gradient le pas utilisé pour le calcul de l'approximation du gradient
% par différences finies.
% max_iter le nombre maximal d'itération de SQP. 
% nfonc_max le nombre d'appel maximal de F.
% borne_inf les bornes inférieures sur la variable, borne_sup ses bornes
% supérieures. Les deux arguments sont optionnels si il n'y a pas 
% de borne sur la variable mais les deux doivent l'être sinon.
% sol le minimiseur obtenu et f_sol son image.
% tableau le tableau des résultats de chaque itération.


%%% Fonctions anonymes

% Gradient du lagrangien
grad_L=@(grad_F_x, lambda) grad_F_x(:,1)+grad_F_x(:,2:end)*lambda; 

% Fonction mérite
merite=@(rho,F_x, norme_c) F_x(1)+rho*norme_c;

% Dérivée partielle de la fonction mérite selon une direction de descente d
grad_merite=@(rho,grad_F_x,norme_c, d) grad_F_x(:,1).'*d-rho*norme_c; 


%%% Initialisation des compteurs, du paramètre de pénalisation et de F(x)

iter=0; 
nfonc=0;
rho=1;
F_x=F(x);
nfonc=nfonc+1;

%%% Initialisations artificielles pour pouvoir exécuter la boucle
%%% au moins une fois.

F_x_old=F_x;
F_x_old(1)=F_x_old(1)+1+tol;
d=zeros(length(x),1);
d(1)=1+tol;
norme_grad_L=1+tol;
grad_L_x_old=0;
H=0;
lambda=0;


%%% Boucle de l'algorithme

while  norme_grad_L>tol && norm(d,1)>tol && norm(F_x-F_x_old,1)> tol && nfonc < max_nfonc && iter<max_iter
    

    %%% Calcul des gradients au point courant
   
    grad_F_x=gradient(F, x, F_x, pas_gradient);
    nfonc=nfonc+length(x);
    grad_L_x=grad_L(grad_F_x, lambda);
    norme_grad_L=norm(grad_L_x, 1);


    %%% Récolte des données du point courant
    
    if iter==0
        tableau=recolte_donnee(iter, nfonc, x, F_x, "None", "None", rho);
    else
        new_line=recolte_donnee(iter,nfonc,x,F_x,lambda,norme_grad_L, rho);
        [tableau]=[tableau; new_line];
    end

    %%% Approximation de la hessienne. 

    if methode=="SR1"
        H=SR1(iter, grad_L_x, grad_L_x_old, H, d);
    elseif methode=="BFGS" 
        H=BFGS(iter, grad_L_x, grad_L_x_old, H, d);
    end
    H=modif_defpos(H);
    

    %%% Résolution du problème quadratique. 

    [new_d, lambda]=KKT_quad_lin(H,F_x, grad_F_x);


    %%% Faire en sorte que d soit une direction de descente de la fonction mérite
    
    norme_c=norm(F_x(2:end), 1);
    check_descente=false;
    check_reset=false;
    has_reset=false;
    while ~check_descente
        grad_merite_x=grad_merite(rho,grad_F_x,norme_c, new_d);
        [H, rho, check_reset, check_descente]=is_descente(grad_merite_x, H, rho, check_reset);
        if check_reset && ~has_reset
            [new_d,lambda]=KKT_quad_lin(H,F_x, grad_F_x);
            has_reset=true;
        end
    end
    d=new_d;
   

    %%% Sauvegarde de F_x et calcul du gradient du Lagrangien au point courant et pour le nouveau
    %%% multiplcateur afin de pouvoir calculer la hessienne à la prochaine
    %%% itération.
    
    F_x_old=F_x;
    grad_L_x_old=grad_L(grad_F_x, lambda); 

    %%% Recherche linéaire

    merite_x=merite(rho,F_x,norme_c);
    if nargin == 7
        [x_pas, nfonc]=armijo(F, merite, rho, d, x, merite_x, grad_merite_x,nfonc);
    else
        [x_pas, nfonc]=armijo(F, merite, rho, d, x, merite_x, grad_merite_x,nfonc, borne_inf,borne_sup);
    end
    d=x_pas-x;
    x=x_pas;


    %%% Fin de l'itération
    
    F_x=F(x);
    nfonc=nfonc+1;
    iter=iter+1;
end


%%% Récolte des données du dernier point
    
grad_F_x=gradient(F, x, F_x, pas_gradient);
nfonc=nfonc+length(x);
grad_L_x=grad_L(grad_F_x, lambda);
norme_grad_L=norm(grad_L_x, 1);
new_line=recolte_donnee(iter,nfonc,x,F_x,lambda,norme_grad_L, rho);
[tableau]=[tableau; new_line];
sol=x;
F_sol=F_x;
end
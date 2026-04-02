function [new_H]=BFGS(iter, grad_L_x, grad_L_x_old, H, d)
% iter l'itération courante.
% grad_L_x la valeur du gradient du lagrangien au point et multiplicateur
% courants.
% grad_L_x_old la valeur du gradient du lagrangien au point précédent et
% au multiplicateur courant.
% H l'appoximation de la hessienne à l'itération k-1.
% d le déplacement du point de départ à l'itération k-1 vers le
% new_H l'appoximation de la hessienne à l'itération k.

if iter==0 
    new_H=eye(length(d));
else
    y=grad_L_x-grad_L_x_old;
    condition=y.'*d;
    if condition > 0 
        terme=H*d;
        new_H=H + y*(y.')/condition - terme*(terme.')/(d.'*terme);
    else 
        new_H=H;
    end
end
end
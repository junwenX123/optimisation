function [x_sol, f_sol, iter] = newton_ralphson(f, x0, tol, max_iter, pas_gradient)
% f la fonction à laquelle on applique les intérations de Newton
% x0 le point de départ. 
% tol la tolérance sur l'erreur de la solution approchée.
% max_iter le nombre d'itérations maximal.
% pas_gradient les pas à utiliser pour les calculs de gradients par
% différences finies.
% x_sol la solution approchée obtenue.
% f_sol son image.
% iter le nombre d'itérations réalisées.

x=x0;
f_x=f(x);
jac_f_x=gradient(f, x, f_x, pas_gradient).';
iter=0;
while tol <norm(f_x,2) && iter<max_iter
    x=linsolve(jac_f_x,-f_x)+x;
    f_x=f(x);
    jac_f_x=gradient(f, x, f_x, pas_gradient).';
    iter=iter+1;
end 
x_sol=x;
f_sol=f_x;
end

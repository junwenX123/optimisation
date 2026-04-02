function [x,lambda]=KKT_quad_lin(H,F_x, grad_F_x)      
% Minimisation d'une fonctionnelle quadratique de forme
% (1/2)<Qx,Qx> + <g,x> sous la contrainte linéaire Ax-b=0 (où A~=0). 
% On suppose que Q est symétrique définie positive et on reformule 
% les conditions KKT sous la contrainte sous la forme d'un système
% linéaire (voir ci-dessous pour le lien avec les arguments).
% x la solution du problème et lambda le vecteur (colonne) des 
% multiplicateurs de Lagrange correspondants.

% Récupération de Q, g, A et b

Q=H;
g=grad_F_x(:,1);
A=grad_F_x(:,2:end).';
b=-F_x(2:end);


% Résolution du problème quadratique linéaire

n=length(g);
m=length(b);
M=[Q,A.';A,zeros(m,m)];
B=[-g;b];
sol=linsolve(M,B);
x=sol(1:n);
lambda=sol(n+1:n+m);
end
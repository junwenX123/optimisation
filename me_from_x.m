function [me, Mtot] = me_from_x(x, k, mu)
% x : rapports Mi/Mf par etage .
% k : indices constructifs par etage .
% mu : masse utile.
% me : masses d'ergols par etage.
% Mtot: masse totale de la fusée.

n = length(x);
Mi_next = mu; % Mi_{j+1}
Mi = zeros(n,1);
Mf = zeros(n,1);
for j = n:-1:1
    denom=(1+k(j))/x(j)-k(j);
    Mi(j) = Mi_next/denom; % Formule établie lors de l'étude analytique
    Mf(j)=Mi(j)/x(j);
    Mi_next = Mi(j);
end
me=Mi-Mf;
Mtot=Mi(1);
end
function [F]=f_ariane1(me, mu, k, ve, delta_requis)
% me le vecteur des masses d'ergol des étages.
% mu la masse utile.
% k le vecteur des indices constructifs des étages.
% ve le vecteur des vitesses d'éjection par étage.
% delta_requis l'incrément de vitesse requis

size=length(me);
M_i=zeros(size+1,1);
M_i(4)=mu;
M_f=zeros(size,1);
M_i(4)=mu;
m_s=k.*me;
delta_v=0;
for i=0:size-1
    M_f(size-i)=M_i(size-i+1)+m_s(size-i);
    M_i(size-i)=M_f(size-i)+me(size-i);
    delta_v=delta_v+ve(size-i)*log(M_i(size-i)/M_f(size-i));
end
F=[M_i(1); delta_v-delta_requis];
end
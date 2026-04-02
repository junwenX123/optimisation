function [F_me]=f_etagement(me, mu, k, ve, vp)
% me le vecteur des masses d'ergol des étages.
% mu la masse utile.
% k le vecteur des indices constructifs des étages.
% ve le vecteur des vitesses d'éjection par étage.
% vp la vitesse propulsive

size=length(me);
M_i=zeros(size+1,1);
M_i(4)=mu;
M_f=zeros(size,1);
m_s=k.*me;
v=0;
for i=0:size-1
    M_f(size-i)=M_i(size-i+1)+m_s(size-i);
    M_i(size-i)=M_f(size-i)+me(size-i);
    v=v+ve(size-i)*log(M_i(size-i)/M_f(size-i));
end
F_me=[mu/M_i(1); v-vp];
end
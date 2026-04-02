function []=tracer_traj(time, result, Mf)
% La fonction trace les graphiques d'intérêt selon les données time,
% result et Mf du simulateur de trajctoire

Rt=6378137;

masse=[result(:, 5);Mf];
altitude=zeros(length(time), 1);
vitesse=zeros(length(time), 1);
for i=1:length(time)
    altitude(i)=(norm(result(i,1:2), 2)-Rt)*10^(-3);
    vitesse(i)=norm(result(i,3:4), 2);
end
abs_altitude=result(:,1);
ord_altitude=result(:,2);

hold on

subplot(2,2,1)
plot([time;time(end)], masse)
xlabel("temps (en s)")
ylabel("masse (en kg)")
title("Evolution de la masse du lanceur")

subplot(2,2,2)
plot(time, altitude)
xlabel("temps (en s)")
ylabel("altitude (en km)")
title("Evolution de l'altitude du lanceur")

subplot(2,2,3)
plot(time, vitesse)
xlabel("temps (en s)")
ylabel("vitesse (en m/s)")
title("Evolution de la vitesse du lanceur")

subplot(2,2,4)
[xunit, yunit]=circle(0,0,6378137,-pi/9, pi/8);
plot(abs_altitude, ord_altitude, '-r', xunit, yunit, '-b')
legend({"fusée", "surface terrestre"},'Location', 'southeast')
title("Trajectoire dans le plan de l'orbite")

hold off

axis equal 
end
function [xunit, yunit] = circle(x,y,r, angle_min, angle_max)
% (x,y) les coordonnées d'un arc de cercle centre d'un cercle de rayon r
% et d'angle angle_max-angle_min.
% xunit les abscisses des points du cercle calculés
% yunit les ordonnées des points du cercle calculés
th = (angle_min:0.02:angle_max).';
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end
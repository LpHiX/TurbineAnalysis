clear; clc; clf; close all
[polars, contours] = readdata( ...
     [ ...
     "sg6043", ...
     "naca23024" ...
     ]);

polar = polars{2};


scatter(polar.Alpha,polar.Cl, 'x');
hold on

%Define t--------------------------
t_index = (-2 < polar.Alpha) & (polar.Alpha < 6);
t_line = polyfit(polar.Alpha(t_index),polar.Cl(t_index),1);
t = @(x) polyval(t_line, x);
plot(-10:10,t(-10:10))

%Define s -------------------------
Cd90 = 2;
Cl90 = 0.08;

b = @(x) x - 56.7 * 0.08 .* sind(x);
s = @(x) (1 + polar.Cl(polar.Alpha==0) * 1.414 * sind(x)) * Cd90 .* sind(b(x)) .* cosd(b(x));
plot(-90:90,s(-90:90))

%Define a1 and a2
a1i = length(polar.Alpha) -1;
diff = abs(polar.Cl - t(polar.Alpha)) > 0.1 & polar.Alpha>0;
scatter(polar.Alpha(diff),polar.Cl(diff), 'o');
a2i = round(0.5 * find(diff, 1, 'first') + 0.5 * a1i);
polar.Alpha(find(diff, 1, 'first'))

scatter([polar.Alpha(a1i), polar.Alpha(a2i)], [polar.Cl(a1i), polar.Cl(a2i)]);

f1 = (polar.Cl(a1i) - s(polar.Alpha(a1i))) / (t(polar.Alpha(a1i)) - s(polar.Alpha(a1i)));
f2 = (polar.Cl(a2i) - s(polar.Alpha(a2i))) / (t(polar.Alpha(a2i)) - s(polar.Alpha(a2i)));

G = ((1 / f1 - 1)/(1 / f2 - 1))^0.25;
a_M = (polar.Alpha(a1i) - G * polar.Alpha(a2i)) / (1 - G);
k = (1 / f2 - 1) * (polar.Alpha(a2i) - a_M).^-4;

f = @(x) 1 ./ (1 + k .* x.^4);
exCl = @(x) f(x) .* t(x) +  (1 - f(x)).* s(x);
plot(-90:90, f(-90:90))
plot(-90:0.1:90, exCl(-90:0.1:90))

grid
legend
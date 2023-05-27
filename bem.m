clear; clc; clf; close all

[polars, contours] = readdata( ...
     [ ...
     "sg6043", ...
     "naca23024" ...
     ]);

%Design Variables --------------------------
B = 3;
R = 0.25;
R0 = 0.025;
N_sections = 20;
design_tsr = 6;

airfoil_profile(1:20) = 1;

r = linspace(R0, R, N_sections);
[d_chord, d_twist] = designblade(B,R,R0,N_sections,design_tsr,airfoil_profile,polars);

N_it = 30;
TSR = 6;
[a1, a2, phi, Cn, Ct] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist);

v_inf = 10;
[Cp,pz,py] = getCp(a1, a2, phi, Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
fprintf("CP is %8.3f, and at TSR %4.2f", Cp, TSR)

TSR_list = 0:0.1:10;
CPlist = zeros(size(TSR_list));
for TSR = TSR_list
    [a1, a2, phi, Cn, Ct] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist);
    [Cp,pz,py] = getCp(a1, a2, phi, Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
    CPlist(TSR==TSR_list) = Cp;
end
plot(TSR_list,CPlist);
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
N_sections = 100;
design_tsr = 5;

airfoil_profile(1:25) = 2;
airfoil_profile(26:100) = 1;

r = linspace(R0, R, N_sections);
[d_chord, d_twist] = designblade(B,R,R0,N_sections,design_tsr,airfoil_profile,polars);

N_it = 30;
TSR = 4.9;
[a1, a2, phi, Cl, Cd, Cn, Ct, aoa, Q] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist);
a1final = a1(end,:);
a2final = a2(end,:);
phifinal = phi(end,:);


v_inf = 10;
[Cp,pz,py,c_ind] = getCp(a1final, a2final, phifinal, Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
fprintf("CP is %8.3f, and at TSR %4.2f \n", Cp, TSR)

v_inf = 10;
torque = startingTorque(airfoil_profile, polars, d_chord,d_twist, B,R,R0,N_sections, v_inf);
fprintf("Starting Torque is %8.3f \n", torque);

TSR_list = 0:0.1:10;
CPlist = zeros(size(TSR_list));
a1list = zeros([length(TSR_list) N_sections]);
for TSR = TSR_list
    [a1, a2, phi, Cl, Cd, Cn, Ct, aoa] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist);
    a1 = a1(end,:);
    a1list(TSR==TSR_list,:) = a1;
    a2 = a2(end,:);
    phi = phi(end,:);
    [Cp,pz,py] = getCp(a1, a2, phi, Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
    CPlist(TSR==TSR_list) = Cp;
end
% plot(TSR_list(CPlist>0),CPlist(CPlist>0))
% plot(0:0.01:10,spline(TSR_list,CPlist,0:0.01:10))
plot(TSR_list,CPlist)
[maxCp, maxtsr] = max(CPlist);
fprintf("CP is %8.3f, and at TSR %4.2f \n", maxCp, TSR_list(maxtsr))
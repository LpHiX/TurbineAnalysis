clear; clc; clf; close all

%1. Example of how to load in airfoil data
airfoilnames = [
     "sg6043"
     "naca23024"
     "sg6040"
     "naca2410"
     "naca4415"
     "naca4421"
     "naca23021"
     "naca632615"
];
[polars, contours] = readdata(airfoilnames);

%2. Example of how to generate a blade geometry (Ben, pls improve)
B = 2;
R = 0.25;
R0 = 0.025;
N_sections = 20;
design_tsr = 5;

%Example of how to change airfoil profile, numbering defined by readdata
airfoil_profile(1:5) = 2;
airfoil_profile(6:20) = 1;

[d_chord, d_twist] = designblade(B,R,R0,N_sections,design_tsr,airfoil_profile,polars);

%3. BEM code to find a1 and a2
N_it = 400;                                 %LOWER THIS TO RUN FASTER
relax_factor = 0.5;                         %High relax, high iterations
TSR = 6; %TSR to test at.

[a1, a2, phi, Cl, Cd, Cn, Ct, aoa, Q] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist,relax_factor);
%BEM code saves all iteration history for debug, only last is useful.
a1final= a1(end,:);
a2final = a2(end,:);
phifinal= phi(end,:);

%4. Use BEM data to calculate Cp, pz and py.
v_inf = 10;
[Cp,pz,py,c_ind] = getCp(a1final, a2final, phifinal, Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
fprintf("CP is %8.3f, and at TSR %4.2f \n", Cp, TSR)

%5. Blade geometry can also be used to find starting torque.
v_inf = 10;
torque = startingTorque(airfoil_profile, polars, d_chord,d_twist, B,R,R0,N_sections, v_inf);
fprintf("Starting Torque is %8.3f \n", torque);

%% That is everything that needs to be used. They can be used to generate plots, for example you can calculate the power curve:

TSR_list = 0:0.1:10;
CPlist = zeros(size(TSR_list));
for TSR = TSR_list
    [a1, a2, phi, Cl, Cd, Cn, Ct, aoa] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist,relax_factor);
    [Cp,pz,py,c_ind] = getCp(a1(end,:), a2(end,:), phi(end,:), Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
    CPlist(TSR==TSR_list) = Cp;
end
plot(TSR_list,CPlist)
[maxCp, maxtsr] = max(CPlist);
fprintf("CP is %8.3f, and at TSR %4.2f \n", maxCp, TSR_list(maxtsr))

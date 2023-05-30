clear; clc; clf; close all

%1
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

for i = 1:length(airfoilnames)
    plot(polars{i}.Alpha, polars{i}.Cl, "DisplayName",airfoilnames(i))
    hold on
end
legend

% %2
% B = 2;
% R = 0.25;
% R0 = 0.025;
% N_sections = 20;
% design_tsr = 5;
% v_inf = 10;
% N_it = 20;
% relax_factor = 0;
% 
% for i = 3
%     airfoil_profile(1:20) = i;
%     [d_chord, d_twist] = designblade(B,R,R0,N_sections,design_tsr,airfoil_profile,polars);
%     TSR_list = 0:0.1:10;
%     CPlist = zeros(size(TSR_list));
%     for TSR = TSR_list
%         %3
%         [a1, a2, phi, Cl, Cd, Cn, Ct, aoa] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist, relax_factor);
%         %4
%         [Cp,pz,py,c_ind] = getCp(a1(end,:), a2(end,:), phi(end,:), Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf);
%         CPlist(TSR==TSR_list) = Cp;
%     end
%     plot(TSR_list(CPlist>0),CPlist(CPlist>0),"DisplayName",airfoilnames(i))
%     hold on
% end
% legend
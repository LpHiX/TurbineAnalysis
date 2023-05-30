B_list = 2:3;
R = 0.25;
R0 = 0.025;
N_sections = 20;
chord = zeros([length(B_list) N_sections]);
for B = B_list
    design_tsr = 5;
    
    airfoil_profile(1:20) = 1;
    
    r = linspace(R0, R, N_sections);
    [d_chord, d_twist] = designblade(B,R,R0,N_sections,design_tsr,airfoil_profile,polars);
    blade(B,:) = d_chord; 
    plot(r, d_chord)
    hold on
end

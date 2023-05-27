function [d_chord, d_twist] = designblade(B,R,R0,N_sections,design_tsr,airfoil_profile,polars)
    r_list = linspace(R0, R, N_sections);
    d_chord = zeros([1,N_sections]);
    d_twist = zeros([1,N_sections]);
    d_aoa = zeros([1,N_sections]);
    for i = 1:N_sections
        r = r_list(i);
        polar = polars{airfoil_profile(i)};
        [~, max_clcd_index] = max(polar.Cl ./ polar.Cd);
        d_aoa(i) = polar.Alpha(max_clcd_index);
        d_cl = polar.Cl(max_clcd_index);
        d_lr = design_tsr * r / R;
        d_phi = atan(2 / (3 * d_lr));
        d_chord(i) = (8 * pi * r * sin(d_phi)) / (3 * B * d_cl * d_lr);
        d_twist(i) = atan(2 / (3 * d_lr)) - deg2rad(d_aoa(i));
    end
end
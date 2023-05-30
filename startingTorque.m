function torque = startingTorque(airfoil_profile, polars, d_chord,d_twist, B,R,R0,N_sections, v_inf)
    r = linspace(R0, R, N_sections);
    aoa = rad2deg(pi/2 - d_twist);
    
    Cl = zeros([1,N_sections]);
    Cd = zeros([1,N_sections]);
    for i = 1:N_sections
        polar = polars{airfoil_profile(i)};
        Cl(i) = interp1(polar.Alpha, polar.Cl, aoa(i), 'linear');
        Cd(i) = interp1(polar.Alpha, polar.Cd, aoa(i), 'linear');
    end
    Cl(isnan(Cl)) = 0;
    density = 1.225;    
    
    p_y = 0.5 * density * v_inf.^2 .* d_chord .* Cl;
    dM = B * p_y .* r;
    torque = trapz(r, dM);
end
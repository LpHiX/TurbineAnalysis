function [a1, a2, phi, Cn, Ct] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist)
    r = linspace(R0, R, N_sections);
    lambda_r = TSR .* (r / R);
    
    %initial angle of attacks
    aoa = zeros([1,N_sections]);

    for i = 1:N_sections
        polar = polars{airfoil_profile(i)};
        [~, max_clcd_index] = max(polar.Cl ./ polar.Cd);
        aoa(i) = polar.Alpha(max_clcd_index);
    end
    
    phi_iterations = zeros([N_it, N_sections]);
    phi_iterations(1, :) = atan(2 ./ (3 * lambda_r));
    a1_iterations = zeros([N_it, N_sections]);
    a2_iterations = zeros([N_it, N_sections]);
    
    sigma_r = (B * d_chord) ./ (2 * pi * r);
    mu = r / R;

    for it = 2:N_it
        %Initialising induction factors
        a1 = a1_iterations(it - 1, :);
        a2 = a2_iterations(it - 1, :);
        phi = phi_iterations(it - 1, :);
        
        %Calculating Cl and Cd for each section
        Cl = zeros([1, N_sections]);
        Cd = zeros([1, N_sections]);
        for i = 1:N_sections
            polar = polars{airfoil_profile(i)};
            Cl(i) = interp1(polar.Alpha, polar.Cl, aoa(i), 'linear');
            Cd(i) = interp1(polar.Alpha, polar.Cd, aoa(i), 'linear');
        end

        %Getting coefficients at better reference frame
        Cn = Cl .* cos(phi) + Cd .* sin(phi);
        Ct = Cl .* sin(phi) - Cd .* cos(phi);
        
        % Iteration starts -------------
        
        Q = 2 / pi * acos(exp(-B ./ (2) * ((1 - mu)./mu) .* sqrt(1 + (lambda_r).^2./(1-a1).^2)));   
        ac = 0.2;
        a1_1 = 1 ./ (1 + (Q .* 4 .* sin(phi).^2) ./ (sigma_r .* Cn));
        K = 4 * sin(phi).^2 .* Q ./ (sigma_r .* Cn);
        a1_2 = 0.5 * (2 + K .* (1 - 2 * ac) - sqrt((K * (1 - 2 * ac) + 2).^2 + 4 * (K * ac^2 - 1)));
        a1(a1<ac) = a1_1(a1<ac);
        a1(a1>ac) = real(a1_2(a1>ac));
    
        a2 = 1 ./ ((Q .* 4 .* cos(phi) .* sin(phi) ./ (sigma_r .* Ct)) - 1);
    
        a1_iterations(it, :) = a1;
        a2_iterations(it, :) = a2;
        phi = atan((1 - a1) ./ ((1 + a2) .* lambda_r));
        phi_iterations(it, :) = phi;
        aoa = rad2deg(phi - d_twist);
    end
end
function [a1, a2, phi, Cl, Cd, Cn, Ct, aoa, Q] = iteratefactors(B,R,R0,N_sections, airfoil_profile,polars, TSR, N_it,d_chord,d_twist)
    r = linspace(R0, R, N_sections);
    lambda_r = TSR .* (r / R);
    
    %initial angle of attacks
    
    phi= zeros([N_it, N_sections]);
    a1 = zeros([N_it, N_sections]);
    a2 = zeros([N_it, N_sections]);
    Cl = zeros([N_it, N_sections]);
    Cd = zeros([N_it, N_sections]);
    Q = zeros([N_it, N_sections]);
    aoa = zeros([N_it,N_sections]);

    phi(1, :) = atan(2 ./ (3 * lambda_r));
    a1(1, :) = 1/3;
    a2(1, :) = 0;

    for i = 1:N_sections
        polar = polars{airfoil_profile(i)};
        [~, max_clcd_index] = max(polar.Cl ./ polar.Cd .* (polar.Alpha < 20 & polar.Alpha > 0));
        aoa(1,i) = polar.Alpha(max_clcd_index);
    end

    sigma_r = (B * d_chord) ./ (2 * pi * r);
    mu = r / R;

    for it = 2:N_it
        for i = 1:N_sections
            polar = polars{airfoil_profile(i)};
            Cl(it,i) = interp1(polar.Alpha, polar.Cl, aoa(it-1,i), 'linear');
            Cd(it,i) = interp1(polar.Alpha, polar.Cd, aoa(it-1,i), 'linear');
        end

        %Getting coefficients at better reference frame
        Cn = Cl(it,:) .* cos(phi(it-1, :)) + Cd(it,:) .* sin(phi(it-1, :));
        Ct = Cl(it,:) .* sin(phi(it-1, :)) - Cd(it,:) .* cos(phi(it-1, :));
        
        % Iteration starts -------------
        
        Q(it,:) = 2 / pi * acos(exp(-B ./ (2) .* ((1 - mu)./mu) .* sqrt(1 + (lambda_r).^2 ./ (1-a1(it-1,:)).^2)));   
        ac = 0.2;
        K = 4 * sin(phi(it-1, :)).^2 .* Q(it,:) ./ (sigma_r .* Cn);
        % a1(it,:) = 1 ./ (K + 1) .* (a1(it-1,:) < ac) + abs(0.5 * (2 + K .* (1 - 2 * ac) - sqrt((K * (1 - 2 * ac) + 2).^2 + 4 * (K * ac^2 - 1)))) .* (a1(it-1,:) >= ac);
        % a2(it,:) = 1 ./ ((Q(it,:) .* 4 .* cos(phi(it-1, :)) .* sin(phi(it-1, :)) ./ (sigma_r .* Ct)) - 1);

        % Relaxation
        a1f = 1 ./ (K + 1) .* (a1(it-1,:) < ac) + abs(0.5 * (2 + K .* (1 - 2 * ac) - sqrt((K * (1 - 2 * ac) + 2).^2 + 4 * (K * ac^2 - 1)))) .* (a1(it-1,:) >= ac);
        a2f = 1 ./ ((Q(it,:) .* 4 .* cos(phi(it-1, :)) .* sin(phi(it-1, :)) ./ (sigma_r .* Ct)) - 1);
        
        relax_factor = 0.5;
        a1(it,:) = relax_factor * a1(it-1,:) + (1 - relax_factor) * a1f;
        a2(it,:) = relax_factor * a2(it-1,:) + (1 - relax_factor) * a2f;

        %--

        phi(it,:) = atan((1 - a1(it,:)) ./ ((1 + a2(it,:)) .* lambda_r));
        aoa(it,:) = rad2deg(phi(it,:) - d_twist);
    end
end
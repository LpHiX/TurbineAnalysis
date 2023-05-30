function [Cp, p_z, p_y, c_ind] = getCp(a1, a2, phi, Cn, Ct, TSR, R, R0, B, N_sections, d_chord, v_inf)
    Cp = 0;
    p_z = 0;
    p_y = 0;

    r = linspace(R0, R, N_sections);
    lambda_r = TSR .* (r / R);
    reference_a2 = -0.5 + 0.5 * sqrt(1 + (4 * a1 .* (1 - a1)) ./ lambda_r.^2);

    error = a2 - reference_a2;
    c_ind = abs(error) < 1e-1;
    % c_ind = ~isnan(a1);
    if length(c_ind) >= 2
        mu = r / R;
        density = 1.225;

        Q = 2 / pi * acos(exp(-B ./ (2) * ((1 - mu)./mu) .* sqrt(1 + (lambda_r).^2./(1-a1).^2)));   
        p_z = 0.5 * density * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Cn .* Q;
        p_y = 0.5 * density * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Ct .* Q;

        M = B * p_y .* r;
        maxP = 0.5 * density * pi * R^2 * v_inf^3;
        Cp = TSR * v_inf / R * trapz(r(c_ind),M(c_ind)) / maxP;
    end
    % mu = r / R;
    % density = 1.225;
    % 
    % Q = 2 / pi * acos(exp(-B ./ (2) * ((1 - mu)./mu) .* sqrt(1 + (lambda_r).^2./(1-a1).^2)));   
    % p_z = 0.5 * density * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Cn .* Q;
    % p_y = 0.5 * density * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Ct .* Q;
    % 
    % M = B * p_y .* r;
    % maxP = 0.5 * density * pi * R^2 * v_inf^3;
    % Cp = TSR * v_inf / R * trapz(r,M) / maxP;
end
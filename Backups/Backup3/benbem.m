clear; clc; clf; close all


polar = getpolars("xf-sg6043-il-50000.csv");


%Initializing Variables --------------------------
B = 3;
N_sections = 26;
R = 0.25;
R0 = 0.025;

TSRList = 0:0.1:10;
CPList = zeros(size(length(TSRList)));
for TSR = TSRList

    N_it = 30;  % converges by 20
    phi_iterations = zeros([N_it, N_sections]);
    a1_iterations = zeros([N_it, N_sections]);
    a2_iterations = zeros([N_it, N_sections]);
    
    %1
    r = linspace(R0, R, N_sections);
    lambda_r = TSR .* (r / R);
    
    phi = atan(2 ./ (3 * lambda_r));
    phi_iterations(1, :) = phi;
    
    %---------------------
    
    [~, max_index] = max(polar.Cl ./ polar.Cd);
    d_aoa = polar.Alpha(max_index);
    d_cl = polar.Cl(max_index);
    d_tsr = 5;
    d_lr = d_tsr * r / R;
    d_phi = atan(2 ./ (3 * d_lr));
    d_chord = (8 * pi .* r .* sin(d_phi)) ./ (3 * B * d_cl * d_lr);
    %plot(r,d_chord);
    d_twist = atan(2 ./ (3 * d_lr)) - deg2rad(d_aoa);
    
    %---------------------
    
    %3
    sigma_r = (B * d_chord) ./ (2 * pi * r);
    aoa = d_aoa;

    for it = 2:N_it
        a1 = a1_iterations(it- 1, :);
        a2 = a2_iterations(it - 1, :);

        % 4)
        Cl = interp1(polar.Alpha, polar.Cl, aoa, 'linear');
        Cd = interp1(polar.Alpha, polar.Cd, aoa, 'linear');
        Cn = Cl .* cos(phi) + Cd .* sin(phi);
        Ct = Cl .* sin(phi) - Cd .* cos(phi);

        mu = r / R;
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

        
    
        % 5)
        phi = atan((1 - a1) ./ ((1 + a2) .* lambda_r));
        phi_iterations(it, :) = phi;
        
        
        % 6)
        aoa = rad2deg(phi - d_twist);
    end
    
    Cl = interp1(polar.Alpha, polar.Cl, aoa, 'linear');
    Cd = interp1(polar.Alpha, polar.Cd, aoa, 'linear');
    
    reference_a2 = -0.5 + 0.5 * sqrt(1 + (4 * a1 .* (1 - a1)) ./ lambda_r.^2);
    
    error = a2 - reference_a2;
    c_ind = find(error < 1e-1);
    if length(c_ind) >= 2
        
        v_inf = 10;
        density = 1.225;
        
        dL = 0.5 * density * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Cn .* Q;

        U = 10;
        dM = 0.5 * density * B * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Ct .* r .* Q;
        maxP = 0.5 * density * pi * R^2 * v_inf^3;
        CPList(TSR==TSRList) = TSR * v_inf / R * trapz(r(c_ind),dM(c_ind)) / maxP;
    end
end
plot(TSRList,CPList);
[maxCp, maxCpIndex] = max(CPList);
fprintf("Max CP is %8.3f, and at TSR %4.2f", maxCp, TSRList(maxCpIndex))

% axis([3 8 0.22 0.4])|


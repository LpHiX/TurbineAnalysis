clear; clc; clf; close all


polar = getpolars("xf-sg6043-il-50000.csv");


%Initializing Variables --------------------------
B = 3;
N_sections = 20;
R = 0.25;
R0 = 0.025;

TSR = 6;

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
d_tsr = 6;
d_lr = d_tsr * r / R;
d_phi = atan(2 ./ (3 * d_lr));
d_chord = (8 * pi .* r .* sin(d_phi)) ./ (3 * B * d_cl * d_lr);
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
    phi = atan((1 - a1) ./ ((1 + a2) .* lambda_r));
    phi_iterations(it, :) = phi;
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
    
    p_z = 0.5 * density * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Q;

    U = 10;

    p_y = 0.5 * density * B * (v_inf .* (1 - a1) ./ sin(phi)).^2 .* d_chord .* Ct .* Q;
    M = p_y .* r;
    maxP = 0.5 * density * pi * R^2 * v_inf^3;
    Cp = TSR * v_inf / R * trapz(r(c_ind),M(c_ind)) / maxP;
end
fprintf("CP is %8.3f, and at TSR %4.2f", Cp, TSR)

airfoil_data = [1.00000     0.00000;
0.95042     0.01245;
0.90089     0.02398;
0.85127     0.03555;
0.80153     0.04693;
0.75163     0.05800;
0.70159     0.06847;
0.65139     0.07809;
0.60105     0.08665;
0.55058     0.09393;
0.50000     0.09974;
0.44932     0.10384;
0.39857     0.10598;
0.34778     0.10587;
0.29700     0.10331;
0.24625     0.09830;
0.19558     0.09066;
0.14504     0.08010;
0.09473     0.06578;
0.06973     0.05667;
0.04492     0.04560;
0.02050     0.03129;
0.00866     0.02159;
0.00418     0.01634;
0.00205     0.01317;
0.00000     0.00000;
0.00795     -0.01017;
0.01082     -0.01214;
0.01634     -0.01517;
0.02950     -0.02013;
0.05508     -0.02664;
0.08027     -0.03123;
0.10527     -0.03476;
0.15496     -0.03972;
0.20442     -0.04290;
0.25375     -0.04460;
0.30300     -0.04499;
0.35222     -0.04407;
0.40143     -0.04172;
0.45068     -0.03814;
0.50000     -0.03356;
0.54942     -0.02823;
0.59895     -0.02239;
0.64861     -0.01629;
0.69841     -0.01015;
0.74837     -0.00430;
0.79847     0.00083;
0.84873     0.00483;
0.89911     0.00704;
0.94958     0.00651;
1.00000     0.00000];

% Step 1: Define the coordinates of the outer edge of the shape
x_airfoil = airfoil_data(:, 1);
y_airfoil = airfoil_data(:, 2);

E = 25e9;
rho = 1240;
omega = TSR * U;

p_z(isnan(p_z)) = 0;
p_y(isnan(p_z)) = 0;

% inertia of airfoil about principle axis along the radius
[I_1,I_2,I_z,I_y,I_zy,alpha,A] = inertia_r(x_airfoil,y_airfoil,d_chord,r);

alpha = deg2rad(alpha); % alpha is B+v

% centrifugal forces
m = rho.*coordinate_integration(r,A,2); % mass in terms of r. the order of the ploynomial is an estimate 
F_centrifugal = rho*omega*coordinate_integration(r,A.*r,2); % centifugal force in terms of r. Check how does the Area*local radius vary to get the degree of polynomial needed  
stress_centrifugal = F_centrifugal./A; % a very rough estimate of centrifugal stresses assuming even for distribution on whole of local airfoil cross section


i_T_z = -p_z+m.*a_z; % y coordinates to be integrated with give x coordinates to find T_z
T_z = coordinate_integration(r,i_T_z,2); % find the value of T_z by integrating the expression:
                                         % -pz(x) + m(x)*az(x) with respect to x

i_T_y = -p_y+m.*a_y; % y coordinates to be integrated with give x coordinates to find T_z
T_y = coordinate_integration(r,i_T_y,2); % find the value of T_z by integrating the expression:
                                         % -pz(x) + m(x)*az(x) with respect to x

M_y = coordinate_integration(r,T_z,3); % finding the bending moments along the blade (radially)
M_z = coordinate_integration(r,-T_y,3); % finding the bending moments along an aerofoil section

% transferring bending moment to principle axis
M_1 = M_y.*cos(alpha)-M_z.*sin(alpha);
M_2 = M_y.*sin(alpha)+M_z.*cos(alpha);

% curvatures about the principle axis from simple beam theory
K_1 = M_1./(E.*I_1);
K_2 = M_2./(E.*I_2);

% transferring the curvatures back to the x and y axis
K_z = -K_1.*sin(alpha)+K_2.*cos(alpha);
K_y = K_1.*cos(alpha)+K_2.*sin(alpha);

% angular deformations
O_y = coordinate_integration(r,K_y,4);
O_z = coordinate_integration(r,K_z,4);

% deflections
u_z = coordinate_integration(r,-O_y,5);
u_y = coordinate_integration(r,O_z,5);

% 
% % bending stresses
% bending_stress_r = bending_stress(M_z,M_y,I_z,I_y,I_zy,x_airfoil,y_airfoil,chord,r);
% 
% 
% subplot(3,2,1);
% plot(r,u_z)
% 
% subplot(3,2,2); 
% plot(r,u_y)
% 
% subplot(3,2,3);
% plot(r,T_z)
% 
% subplot(3,2,4); 
% plot(r,T_y)
% 
% subplot(3,2,5);
% plot(r,stress_centrifugal)
% 
% subplot(3,2,6); 
% plot(r,bending_stress_r)


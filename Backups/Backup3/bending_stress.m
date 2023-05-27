function bending_stress = bending_stress(M_z,M_y,I_z,I_y,I_zy,x_airfoil_coordinates,y_airfoil_coordinates,chord,r)
% a function to calculate the maximum stress due to bending in a cross
% section of airfoil in terms of r
% inputs: x_airfoil_coordinates - x coordinates of airfoil used
%         x_airfoil_coordinates - y coordinates of airfoil used
%              chord            - chord distribution
%                r              - points in radial direction where the forces and the
%                                 acceleration, etc are evaluated
%           I_1 - moment of inertia about the principle x axis in the
%                 radial direction
%           I_2 - moment of inertia about the principle y axis in the
%                 radial direction
%           M_1 - bending moment about the principle x axis in the radial
%                 direction
%           M_2 - bending moment about the principle y axis in the radial
%                 direction
%       I_z,I_y,I_zy - % moment of inertia about the x' and y' axis which is parellel to the ploting axis and has centre at the centroid of the airfoil

% for this calculations the slides change x to z

% preallocating the matrix
c = (r);
bending_stress = zeros(2,c);

% precalcualtions for efficiency
A = (M_y.*I_z-M_z.*I_zy)./(I_z.*I_y-I_zy.^2);
B = (M_z.*I_y-M_y.*I_zy)./(I_z.*I_y-I_zy.^2);


for i = 1:length(r)
    z_airfoil_section = x_airfoil_coordinates.*chord(i); % calculating new x coordinates of airfoil that relate to the chord at the specific r value
    y_airfoil_section = y_airfoil_coordinates.*chord(i); % calculating new y coordinates of airfoil that relate to the chord at the specific r value
    [z_mesh,y_mesh] = mesh_coordinates(z_airfoil_section,y_airfoil_section); % getting the coordinates of the mesh that lie within the airfoil
    stress_x = -A(i)*z_mesh-B*y_mesh; % calculating stress at each coordinate of the mesh
    positive_max_stress_x = max(stress_x); % finding the max tension (positive stress)
    negative_max_stress_x = min(stress_x); % finding the max compression (negative stress)
    % abs_stress_x = abs(stress_x);
    % max_abs_stress = max(abs_stress_x);
    % [r,c] = find(abs_stress_x == max_abs_stress);
    % negative_max_stress_x = stress_x(r,c); % finding the max compression (negative stress)
    % if negative_max_stress_x > 0
    %     negative_max_stress_x = 0;
    % end
    bending_stress(1,i) = negative_max_stress_x; % max compressive forces
    bending_stress(2,i) = positive_max_stress_x; % max tension forces
end

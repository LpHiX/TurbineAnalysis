function [u_z,u_y,T_z,T_y,stress_centrifugal,bending_stress_r] = structural_deflections(p_z,p_y,a_z,a_y,r,x_airfoil,y_airfoil,E,omega,rho,chord)
% This is a function to calculate deflections of the blade and the shear
% forces acting on the blade as x (radial direction) varies
% inputs:  p_z - varying axial loads in terms of r
%          p_y - varying tangential loads in terms of r
%          a_z - acceleration in axial plane
%          a_y - acceleration in tangential plane
%           r  - points in radial direction where the forces and the
%                acceleration is evaluated
%           E  - young's modulus
%         y_airfoil - y coordinates of an airfoil
%         x_airfoil - x coordinates of an airfoil
%          omega - turbine angular velocity
%           rho - density of the material
%          chord - chord distribution in terms of r
% outputs: T_z - shear stress in axial direction in terms of r
%          T_y - shear stress in tangential direction in terms of r
%          u_z - deflection in axial direction in terms of r
%          u_z - deflection in tangential direction in terms of r
% all the inputs and outputs are column vectors

% inertia of airfoil about principle axis along the radius
[I_1,I_2,I_z,I_y,I_zy,alpha,A] = inertia_r(x_airfoil,y_airfoil,chord,r);

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


% bending stresses
bending_stress_r = bending_stress(M_z,M_y,I_z,I_y,I_zy,x_airfoil,y_airfoil,chord,r);


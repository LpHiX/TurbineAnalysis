function [I_1,I_2,I_x,I_y,I_xy,alpha,A] = inertia_r(x_airfoil_coordinates,y_airfoil_coordinates,chord,r)
% a function to evaluate the the inertia of airfoil sections in the radial
% direction by connsdidering chord distribution. Also the function gives
% the angle between the ploting axis and the principle axis along the radial direction
% inputs: x_airfoil_coordinates - x coordinates of airfoil used
%         x_airfoil_coordinates - y coordinates of airfoil used
%              chord            - chord distribution
%                r              - points in radial direction where the forces and the
%                                 acceleration, etc are evaluated
% output: alpha - angle between the x' axis and principle x axis in the radial direction. Given in
%                 degress
%           I_1 - moment of inertia about the principle x axis in the
%                 radial direction
%           I_2 - moment of inertia about the principle y axis in the
%                 radial direction
%            A  - area of airfoil section in terms of r
%       I_x,I_y,I_xy - % moment of inertia about the x' and y' axis which is parellel to the ploting axis and has centre at the centroid of the airfoil

% preallocating the matrix
sz = size(r);
aplha = zeros(sz);
I_1 = zeros(sz);
I_2 = zeros(sz);
I_x = zeros(sz);
I_y = zeros(sz);
I_xy = zeros(sz);
A = zeros(sz);

for i = 1:length(r)
    x_airfoil_section = x_airfoil_coordinates.*chord(i); % calculating new x coordinates of airfoil that relate to the chord at the specific r value
    y_airfoil_section = y_airfoil_coordinates.*chord(i); % calculating new y coordinates of airfoil that relate to the chord at the specific r value
    airfoil = polyshape(x_airfoil_section,y_airfoil_section);
    [X_E,Y_E] = centroid(airfoil); % calculating the coordinates of centroid of the airfoil about the ploting axis
    A(i) = polyarea(x_airfoil_section,y_airfoil_section); % calculating the area of the airfoil section
    [I_xR,I_yR,I_xyR] = inertia_of_airfoil_section(x_airfoil_section,y_airfoil_section); % calculating the inertia of that section of airfoil about the ploting axis
    I_x(i) = I_xR-(Y_E^2)*A(i); % moment of inertia about the x' and y' axis which is parrellel to the ploting axis and has centre at the centroid of the airfoil
    I_y(i) = I_yR-(X_E^2)*A(i);
    I_xy(i) = I_xyR-X_E*Y_E*A(i);
    alpha = 0.5*atan(2*I_xy(i)/(I_y(i)-I_x(i))); % angle between the x' axis and principle x axis
    I_1(i) = I_x(i)-I_xy(i)*tan(alpha); % moment of inertia about the principle axis
    I_2(i) = I_y(i)-I_xy(i)*tan(alpha);
    alpha(i) = rad2deg(alpha); % converting alpha to degress
end

function [I_xR,I_yR,I_xyR] = inertia_of_airfoil_section(x_airfoil_coordinates,y_airfoil_coordinates)
% calculating inertia for each airfoil cross section
% inputs: x_airfoil_coordinates - x coordinates of airfoil used
%         x_airfoil_coordinates - y coordinates of airfoil used
% outputs: I_xR - moment of inertia of the 2D airfoil about the x axis of 
%                 the airfoil plot
%          I_yR - moment of inertia of the 2D airfoil about the y axis of 
%                 the airfoil plot
%          I_xyR - cross product of moment of inertia of the 2D airfoil about origin of 
%                  the airfoil plot



% Step 1: Define the coordinates of the outer edge of the shape
x_contour =x_airfoil_coordinates;
y_contour = y_airfoil_coordinates;


% Step 3: Determine the bounding box
xmin = min(x_contour);
xmax = max(x_contour);
ymin = min(y_contour);
ymax = max(y_contour);

% Step 4: Create a grid of points within the bounding box
xresolution = 100;
yresolution = 20;
[X, Y] = meshgrid(linspace(xmin, xmax, xresolution), linspace(ymin, ymax, yresolution));



elements = repmat(polyshape, 1, xresolution*yresolution);
x_centroids = zeros(size(X));
y_centroids = zeros(size(X));

for i = 1:xresolution-1
    for j = 1:yresolution-1
        x1 = X(j, i);
        y1 = Y(j, i);
        x2 = X(j, i+1);
        y2 = Y(j, i+1);
        x3 = X(j+1, i+1);
        y3 = Y(j+1, i+1);
        x4 = X(j+1, i);
        y4 = Y(j+1, i);
        x_element = [x1,x2,x3,x4];
        y_element = [y1,y2,y3,y4];
        if inpolygon(x_element, y_element, x_contour, y_contour)
            elements(i*yresolution + j)  = polyshape(x_element,y_element);
            [x_centroid,y_centroid] = centroid(elements(i*yresolution + j));
            x_centroids(i*yresolution + j) = x_centroid; % x coordinate of centroid of each element
            y_centroids(i*yresolution + j) = y_centroid; % y coordinate of centroid of each element
            dA(i*yresolution + j) = polyarea(x_element,y_element); % area of each element
        end
    end
end
A = area(elements); 
elements = elements(A>0);
x_centroids = x_centroids(A>0);
y_centroids = y_centroids(A>0);
dA = dA(A>0);

% calculating the inertia of the cross section of airfoil using basic principles

% trapz might not work due to dA being same hence when plotted their will
% be a vertical line giving the inetegral to be zero
% I_xR = trapz(dA,y_centroids.^2);
% I_yR = trapz(dA,x_centroids.^2);
% I_xyR = trapz(dA,x_centroids.*y_centroids);

% using summation instead of integral
I_xR = sum(dA.*(y_centroids.^2));
I_yR = sum(dA.*(x_centroids.^2));
I_xyR = sum(dA.*x_centroids.*y_centroids);

end
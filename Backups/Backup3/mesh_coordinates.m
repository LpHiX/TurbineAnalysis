function [x_mesh,y_mesh] = mesh_coordinates(x_airfoil_coordinates,y_airfoil_coordinates)
% a function to find the x and y coordinates of the edges of elements of
% mesh of an airfoil. Remember to have the resolution same in this and
% 'inertia_of_airfoil_section' function.
% inputs: x_airfoil_coordinates - x coordinates of airfoil used
%         x_airfoil_coordinates - y coordinates of airfoil used
% ouputs: x_mesh - x coordinate of mesh within the airfoil
%         y_mesh - y coordinates of mesh within the airfoil

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
            
        end
    end
end
A = area(elements); 
elements = elements(A>0);
x_mesh = X(A>0);
y_mesh = Y(A>0);

end

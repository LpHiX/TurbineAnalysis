clear; clc; clf; close all
data = [1.00000     0.00000;
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
x = data(:, 1);
y = data(:, 2);


% Step 2: Create a polyshape object for the shape
shape = polyshape(x, y);

% Step 3: Determine the bounding box
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% Step 4: Create a grid of points within the bounding box
xresolution = 50;
yresolution = 10;
[X, Y] = meshgrid(linspace(xmin, xmax, xresolution), linspace(ymin, ymax, yresolution));

% Step 5: Use inpolygon to check which mesh points lie within the shape
mask = inpolygon(X(:), Y(:), x, y);

% Step 6: Get the coordinates of the mesh points within the shape
mesh_x_within_shape = X(mask);
mesh_y_within_shape = Y(mask);

% Step 7: Plot the shape and the mesh points within the shape
figure;
plot(shape);
hold on;
plot(mesh_x_within_shape(1:10), mesh_y_within_shape(1:10), 'b.', 'MarkerSize', 10);


% Step 8: Set the aspect ratio to equal and adjust the axes
axis equal;

% Step 9: Add labels and title
xlabel('X');
ylabel('Y');
title('Mesh Points within Shape');










%-------------------------------------------------------------------


% mesh grid
X = mesh_x_within_shape;
Y = mesh_y_within_shape;

% Get size of the mesh grid
[nrows, ncols] = size(X);


% Iterate over each element in the mesh grid
for i = 1:nrows-1
    for j = 1:ncols-1
        % Get the coordinates of the four corners of the current element
        x1 = X(i, j);
        y1 = Y(i, j);
        x2 = X(i+1, j);
        y2 = Y(i+1, j);
        x3 = X(i+1, j+1);
        y3 = Y(i+1, j+1);
        x4 = X(i, j+1);
        y4 = Y(i, j+1);
        
        x_element = [x1,x2,x3,x4];
        y_element = [y1,y2,y3,y4];

        element  = polyshape(x,y);
        
        % Add the coordinates of the edges of the current element to the arrays
        [x_centroid,y_centroid] = centroid(element);
        x_centroid(i,j)=x_centroid;
        y_centroid(i,j)=y_centroid;

        % Area of each element
        dA(i,j) = polyarea(x_element,y_element);
    end
end

% Plot the mesh grid and the edges of each element
plot(x_centroid,y_centroid,'g*')

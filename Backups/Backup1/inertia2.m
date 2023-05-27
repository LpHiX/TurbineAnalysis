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
x_contour = data(:, 1);
y_contour = data(:, 2);


% Step 2: Create a polyshape object for the shape
shape = polyshape(x_contour, y_contour);

% Step 3: Determine the bounding box
xmin = min(x_contour);
xmax = max(x_contour);
ymin = min(y_contour);
ymax = max(y_contour);

% Step 4: Create a grid of points within the bounding box
xresolution = 100;
yresolution = 20;
[X, Y] = meshgrid(linspace(xmin, xmax, xresolution), linspace(ymin, ymax, yresolution));


% Step 5: Plot the shape and the mesh points within the shape
figure;
plot(shape);
hold on;
plot(X, Y, 'b.', 'MarkerSize', 10);

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
            x_centroids(i*yresolution + j) = x_centroid;
            y_centroids(i*yresolution + j) = y_centroid;
            fill(x_element,y_element, "r")
        end
    end
end
A = area(elements);
elements = elements(A>0);
x_centroids = x_centroids(A>0);
y_centroids = y_centroids(A>0);
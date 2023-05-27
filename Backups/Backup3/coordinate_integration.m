function I_y = coordinate_integration(x,y,n)
% A function which evaluates the integral of the function represented
% by the x and y coordinates for a give n 
% inputs: x - x coordintaes
%         y - y cocordinates
%         n - the degree of polynomial estimated for the plot of x and y
%             coordinates
% outputs: I_y - the coordinates of the integrated function evaluated at
%                the original (input) x values
% all the inputs and outputs are column vectors. n is a positive scalar more than 1

% Fit a polynomial of degree n to the data points
p = polyfit(x, y, n);

% Create a symbolic variable and function
syms x_sym;
f_sym = sym(poly2sym(p, x_sym));

% Integrate the symbolic function
integral_f_sym = int(f_sym, x_sym);

% Evaluating the integrated function at the original (input) x values
I_y = subs(integral_f_sym,x_sym,x);

end
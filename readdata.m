function [polars, contours] = readdata(folders)
polars = cell(1,length(folders));
contours = cell(1,length(folders));

for i = 1:length(folders)
    %Opens polars
    folder = folders(i);
    polars{i} = load(strcat("airfoils/",folder, "/", folder, ".mat")).finaldata;
    %Opens contours
    contourpath = strcat("airfoils/",folder,"/",folder, ".dat");
    contours{i} = readtable(contourpath, 'HeaderLines', 1);
end
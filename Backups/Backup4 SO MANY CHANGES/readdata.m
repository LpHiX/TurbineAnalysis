function [polars, contours] = readdata(folders)
polars = cell(1,length(folders));
contours = cell(1,length(folders));

for i = 1:length(folders)
    folder = folders(i);

    %Opens polars
    % polarpath = strcat(folder,"/xf-",folder,"-il-100000.csv");
    % 
    % fid = fopen(polarpath, 'rt');
    % lineCount = 0;
    % while ~feof(fid)
    %     line = fgetl(fid);
    %     lineCount = lineCount + 1;
    %     if contains(line, 'Alpha')
    %         break;
    %     end
    % end
    % fclose(fid);
    % polars{i} = readtable(polarpath, 'HeaderLines', lineCount - 1);
    polars{i} = load(strcat(folder, "/data.mat")).finaldata;
    


    %Opens contours
    contourpath = strcat(folder,"/",folder);
    contours{i} = readtable(contourpath, 'HeaderLines', 1);
end
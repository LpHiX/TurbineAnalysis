function polars = getpolars(fname)

fid = fopen(fname, 'rt');
lineCount = 0;
while ~feof(fid)
    line = fgetl(fid);
    lineCount = lineCount + 1;
    if contains(line, 'Alpha')
        break;
    end
end
fclose(fid);
polars = readtable(fname, 'HeaderLines', lineCount - 1);
end
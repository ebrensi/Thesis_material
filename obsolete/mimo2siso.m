function mimio2siso(data)
% mimio2siso(data)
%
% convert MIMO data file to multiple SISO files with identical left
% and right side multipliers.

[A E B] = inputdata(data);
m = size(B,2);

for i = 1:m
    b = B(:,i);
    c = b;
    fname = sprintf('%ss_%d',data,i);
    fprintf('%s\n',fname);
    save(fname,'A','E','c','b');
end


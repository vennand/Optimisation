function [data] = adjust_number_of_interval(data)
for i=1:data.realNint
    Nint(i)= numel(data.frames(1):i:data.frames(end));
end
[~, indexOfMin] = min(abs(Nint-data.Nint));

if data.Nint == Nint(indexOfMin)-1
%     fprintf('The number of interval is: %d\n', data.Nint);
else
    data.Nint = Nint(indexOfMin)-1;
    fprintf('The number of interval has been modified to: %d\n', data.Nint);
end
data.step = indexOfMin;
end
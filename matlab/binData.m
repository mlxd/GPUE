function [g6B,bin] = binData(g6C,binMax,bins)
%Bin the data returned from g6_struct into equi-separated partitions
% g6C:       Matrix containing the g6 values output from g6_struct
% binMax:    The maximum binning value
% bins:      The number of bins to take
%Returns
% g6B:       The binned g6 data
% bin:       The bin values of g6B

bin = linspace(0,binMax,bins);
g6B = zeros(length(bin)-1,1);
for kk=1:(length(bin)-1)
    idx = find( g6C(1,:) >= bin(kk) & g6C(1,:) < bin(kk+1) );
    if length(idx) >0
        g6B(kk) = sum(g6C(2,idx))./length(idx); %Average the values in a bin
    else
        g6B(kk) = 0;
    end
end

end
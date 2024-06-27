function [f, activeInd] = calculateFunctions0(x,th,N)
p = mod(x(1),2*pi);
thc = x(2);

activeInd = find((th > p-thc + 1e-15) & (th < p + thc - 1e-15));

if p + thc >= 2*pi
    ind2 = find(th < p+thc-2*pi);
    activeInd = [activeInd; ind2];
end

if p-thc < 0
    ind2 = find(th > p-thc+2*pi);
    activeInd = [activeInd; ind2];
end

if isempty(activeInd)
    fH0 = -1e3;
    frho = -1e3;
    fpsi = 1e3;
else
    diff = cos(th(activeInd) - p) - cos(thc);
    fH0 = 1/N*sum(diff);
    frho = 1/N*sum(diff.*cos(th(activeInd) - p));
    fpsi = 1/N*sum(diff.*sin(th(activeInd) - p));
end

f = [fH0; frho; fpsi];

end
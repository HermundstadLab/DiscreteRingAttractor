function [f, activeInd] = calculateFunctions(x,N)

% x(:,1) = psi values
% x(:,2) = thc values
% N = number of computational units

th = linspace(0,2*pi,N+1); th(end) = []; th = th';

f = zeros(size(x,1),3);
activeInd = zeros(size(x,1),N);

for i = 1:size(x,1)
    p = mod(x(i,1),2*pi);
    thc = x(i,2);
    
    actInd = find((th > p-thc) & (th < p + thc));
    
    if p + thc >= 2*pi
        ind2 = find(th < p+thc-2*pi);
        actInd = [actInd; ind2];
    end
    
    if p-thc < 0
        ind2 = find(th > p-thc+2*pi);
        actInd = [actInd; ind2];
    end
    
    if isempty(actInd)
        fH0 = NaN;
        frho = NaN;
        fpsi = NaN;
    else
        diff = cos(th(actInd) - p) - cos(thc);
        fH0 = 1/N*sum(diff);
        frho = 1/N*sum(diff.*cos(th(actInd) - p));
        fpsi = 1/N*sum(diff.*sin(th(actInd) - p));
    end
    
    f(i,:) = [fH0; frho; fpsi];
    activeInd(i,actInd) = 1;
    
end

end
function J0 = computeJI(J1Vec,A,C0,psiSamp,N,numSamps)

if nargin < 6
    numSamps = 100;
end

J0 = zeros(length(J1Vec),1);

for j1 = 1:length(J1Vec)
    
    J1 = J1Vec(j1);
    
    if length(psiSamp) <= 2*N
        % find thetac corresponding to psiFP with JE = J1
        thcinit = pi/2;
        thc = zeros(size(psiSamp));
        for i = 1:length(psiSamp)
            [thc(i),fval,flag,~] = fzero(@(thc)([0 J1 0]*calculateFunctions([psiSamp(i) thc],N)' - 1),thcinit);
            if flag < 0
                disp(['Something went wrong, flag for root finder: ' num2str(flag)])
            end
        end
        
        psi = psiSamp;
        
    else
        thcSamp = linspace(pi/N, (N-1)*pi/N, numSamps+1);
        thcSamp(end) = []; thcSamp = thcSamp';
        
        frho = zeros(length(thcSamp),length(psiSamp));
        for p = 1:length(psiSamp)
            [f,~] = calculateFunctions([psiSamp(p)*ones(size(thcSamp)),thcSamp],N);
            frho(:,p) = f(:,2);
        end
        
        Mc = contourc(psiSamp, thcSamp, frho, [1/J1 1/J1]);
        psi = Mc(1,2:end)';
        thc = Mc(2,2:end)';
        
        [f,Nact] = calculateFunctions([psi thc],N);
        J0ub = -cos(thc)./f(:,1);
        
    end
    
    [f,~] = calculateFunctions([psi thc],N);
    f0 = f(:,1);
    % J0 to approx achieve desired amplitude
    J0temp = ((C0/A-1)*cos(thc) - C0/A)./f0;
    J0(j1) = min(J0temp);
    
    % theoretical upper bound on J0
    J0ub = min(-cos(thc)./f(:,1));
    
    if J0ub - J0(j1) <= 0
        disp(['calculated J0 outside of marginally stable regime, JE = ' num2str(J1)])
    end
    
end

end


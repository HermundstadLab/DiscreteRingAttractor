function E = energyLandscape(J1,J0,C0,N,psiSamples,thcSamples)

%% PREFERRED FEATURE SPACE
th = linspace(0,2*pi,N+1); th(end) = [];

%% CONNECTIVITY
W = 1/N*(J0 + J1*cos(th - th'));

%% CALCULATE ENERGY
Erec = NaN(length(thcSamples),length(psiSamples));
Eext = Erec;
Eleak = Erec;

frho = zeros(length(thcSamples), length(psiSamples));
rho = frho;

for p = 1:length(psiSamples)
    [f,~] = calculateFunctions([psiSamples(p)*ones(size(thcSamples)), thcSamples], N);
    f0 = f(:,1);
    frho(:,p) = f(:,2);
    rho(:,p) = -C0./(2*(cos(thcSamples) + J0*f0));
    %validInd = (rho(:,p) > 0);
    %validInd = ones(size(rho(:,p)));
    %validThc = thcSamples(validInd);
    %ErecTemp = zeros(length(validThc),1);
    %EextTemp = ErecTemp;
    %EleakTemp = ErecTemp;
    for c = 1:length(thcSamples)
        h = 2*rho(c,p)*(cos(th-psiSamples(p)) - cos(thcSamples(c)));
        r = poslin(h);
        
%         ErecTemp(c) = -1/2*r*W*r';
%         EextTemp(c) = -C0*sum(r);
%         EleakTemp(c) = 1/2*sum(r.^2);
        
        Erec(c,p) = -1/2*r*W*r';
        Eext(c,p) = -C0*sum(r);
        Eleak(c,p) = 1/2*sum(r.^2); 
    end
    
%     Erec(validInd,p) = ErecTemp;
%     Eext(validInd,p) = EextTemp;
%     Eleak(validInd,p) = EleakTemp;
%     
end

E = Erec + Eext + Eleak;

end


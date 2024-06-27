function dh = dynSys(t,h,W,Wcw,Wccw,C0,tau,v,dt,addNoise)

ind = t/dt + 1;
indf = floor(ind);
indc = ceil(ind);

if length(v) > 1
    if indf == indc
        v = v(ind);
    else
    v = v(indf)*(indc - ind) + v(indc)*(ind - indf);
    end
end

if size(addNoise,2) > 1
    if indf == indc
        addNoise = addNoise(:,ind);
    else
    addNoise = addNoise(:,indf)*(indc - ind) + addNoise(:,indc)*(ind - indf);
    end
end

vcw = poslin(v);
vccw = poslin(-v);

I = C0;

dh = W*poslin(h) + vcw*Wcw*poslin(h) + vccw*Wccw*poslin(h) + addNoise;
dh = dh - h;
dh = dh + I;
dh = dh/tau;

end
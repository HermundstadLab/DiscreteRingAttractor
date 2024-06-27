function [h,H0,rho,psi,thc] = simM(W,Wcw,Wccw,C0,tau,v,IC,t,N,addNoise)

dt = (t(end) - t(1))/(length(t)-1);

[t,h] = ode45(@(t,h)dynSys(t,h,W,Wcw,Wccw,C0,tau,v,dt,addNoise),t,IC);
DFT = fft(h');

H0 = 1/N*DFT(1,:); H0 = H0';
rho = 1/N*abs(DFT(2,:)); rho = rho';
psi = -angle(DFT(2,:)); psi = psi';
thc = acos(-H0./(2*rho));

end


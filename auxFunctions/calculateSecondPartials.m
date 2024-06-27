function H = calculateSecondPartials(rho,thc,psi, J0, J1, C0, th, N)

[f,kact] = calculateFunctions0([psi thc], th,N);
f0 = f(1);
frho = f(2);
fpsi = f(3);

Nact = length(kact);

W = 1/N*(J1*cos(th(kact) - th(kact)') + J0);

s1 = sum((cos(th(kact) - psi) - cos(thc)).^2);
s2 = (cos(th(kact) - psi) - cos(thc))'*W*(cos(th(kact) - psi) - cos(thc));
drho2 = 4*(s1 - s2);

s1 = sum((cos(th(kact) - psi) - cos(thc)));
s2 = sum(W*(cos(th(kact) - psi) - cos(thc)));
drhodthc = 8*rho*sin(thc)*(s1 - s2) - 2*C0*Nact*sin(thc);

s1 = sum(sin(th(kact)-psi));
s2 = N*fpsi;
s3 = sum(W*(cos(th(kact) - psi) - cos(thc)).*sin(th(kact) - psi));
drhodpsi = 2*(-C0*s1 + 4*rho*s2 - 4*rho*s3);

s1 = N*f0;
s2 = sum(W*(cos(th(kact) - psi) - cos(thc)));
s3 = sum(sum(W));
dthc2 = 4*rho^2*cos(thc)*(s1 - s2) + 4*rho^2*(sin(thc))^2*(Nact - s3) -2*C0*Nact*rho*cos(thc);

s1 = sum(sin(th(kact)-psi));
s2 = sum(W*sin(th(kact) - psi));
dthcdpsi = 4*rho^2*sin(thc)*(s1 - s2);

s1 = sum((sin(th(kact) - psi)).^2);
s2 = frho*N;
s3 = (sin(th(kact) - psi))'*W*sin(th(kact) - psi);
s4 = (cos(th(kact) - psi))'*W*(cos(th(kact) - psi) - cos(thc));
s5 = sum(cos(th(kact) - psi));
dpsi2 = 4*rho^2*(s1 - s2 - s3 + s4) + 2*C0*rho*s5;


H = [drho2 drhodthc drhodpsi; drhodthc dthc2 dthcdpsi; drhodpsi dthcdpsi dpsi2];
end
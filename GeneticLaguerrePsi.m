%% Genetic Laguerre Psi
% Mode, initial amplitude and omega_0
n = 0;
m = 0;
A0 = 1;
w0 = 1;
[theta_0,~] = cart2pol(xx_0,yy_0);
[theta_1,~] = cart2pol(xx_1,yy_1);

% Psi base
fprintf('Calculating initial PSI state... ');

rw0_0 = (r_0/w0);
rw0_1 = (r_1/w0);

psi_0 = A0 * (rw0_0).^m .* laguerreL(n,2*rw0_0.^2).^m .* exp(-rw0_0.^2 + 1i*(m*theta_0));
psi_1 = A0 * (rw0_1).^m .* laguerreL(n,2*rw0_1.^2).^m .* exp(-rw0_1.^2 + 1i*(m*theta_1));

psi_base =  psi_0 + psi_1;

fprintf('Done!\n');
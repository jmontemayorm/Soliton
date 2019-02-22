%% Non-Lineal Schrödinger Equation (NLSE)
%
%   dP   1 d^2P   1 d^2P
%  i-- + - ---- + - ---- + |P|^2 P + V P= 0
%   dz   2 dx^2   2 dy^2
%

%% Settings

% Genetic settings
finalGeneration = 3;
populationSize = 6;
polynomialTerms = 5;

% Parallel settings
maxWorkers = 6;

% Space definitions
numOfPoints = 256;

windowSize = 10;
x = linspace(-windowSize,windowSize,numOfPoints);
dx = x(2)-x(1);

[xx,yy] = meshgrid(x,x); % Base grid

[xx_0,yy_0] = meshgrid(x-2,x);
r_0 = sqrt(xx_0.^2+yy_0.^2);

[xx_1,yy_1] = meshgrid(x+2,x);
r_1 = sqrt(xx_1.^2+yy_1.^2);

limitZ = 4;
dz = dx^2/8;
doubleSteps = ceil(length(0:dz:limitZ)/2);

% Mode, initial amplitude and omega_0
n = 0;
m = 0;
A0 = 1;
w0 = 1;
[theta_0,~] = cart2pol(xx_0,yy_0);
[theta_1,~] = cart2pol(xx_1,yy_1);

% Psi base
rw0_0 = (r_0/w0);
rw0_1 = (r_1/w0);
psi_0 = A0 * (rw0_0).^m .* laguerreL(n,2*rw0_0.^2).^m .* exp(-rw0_0.^2 + 1i*(m*theta_0));
psi_1 = A0 * (rw0_1).^m .* laguerreL(n,2*rw0_1.^2).^m .* exp(-rw0_1.^2 + 1i*(m*theta_1));
psi_base =  psi_0 + psi_1;

% Initial potential guess
V0 = -(abs(psi_base).^2) + 4*n + 2*m + 2 - 2*r_0.^2 - 2*r_1.^2;

%% Initial population
V_terms = rand(populationSize,polynomialTerms);

% Pre-calculated polynomials
polyTerms = zeros([size(xx,1),size(xx,2),polynomialTerms]);
for p = 1:polynomialTerms
    polyTerms(:,:,p) = xx.^(p-1);
end

% Initial potentials
V = cell(1,populationSize);
for specimen = 1:populationSize
    tempV = V0;
    for p = 1:polynomialTerms
        tempV = tempV + V_terms(specimen,p) * polyTerms(:,:,p);
    end
    V{specimen} = tempV;
end


%% Evolution
tic
psi_results = cell(1,populationSize);
for generation = 1:finalGeneration
    % Parallel propagation
    parfor (specimen = 1:populationSize,maxWorkers)
        psi_results{specimen} = SolitonPropagator2D(psi_base,V{specimen},windowSize,dz,doubleSteps)
    end
    
    % Evaluation
    
    % Evolution mechanisms
end
toc

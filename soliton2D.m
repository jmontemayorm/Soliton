%% Non-Lineal Schrödinger Equation (NLSE)
%
%   dP   1 d^2P   1 d^2P
%  i-- + - ---- + - ---- + |P|^2 P = 0
%   dz   2 dx^2   2 dy^2
%

%% Start

% Time script
timeScript = 0;

if timeScript == 1
    tic
end

% 0 = exp(-r.^2)
% 1 = N1*sech(r)
% 2 = N1*sech(r1)+N2*sech(r2) % r1 = sqrt((x-x01).^2+(y-y01).^2)
% 3 = sech(r+a).*exp(1i*alpha*r)+sech(r-a).*exp(-1i*alpha*r)
% 4 = sech(r+a) + exp(-(r-a).^2)
% 5 = sum(sech(r@(x0,y0))*exp(2i*pi*alphaVec), momentum p = -1,0,1), circle r0
% 6 = 5, but ellipse with semi-x and semi-y axes sX and sY
% 7 = special case
beamProfile = 5;
numOfSol = 6;
r0 = 3;
alpha = 1;
alphaVec = alpha*ones(1,numOfSol); %linspace(0,2*pi*(1-1/numOfSol),numOfSol);
alphaVec(2:2:end) = 0; % Alternate values
sX = 5; % Horizontal "radius" (ellipse)
sY = 4; % Vertical "radius" (ellipse)
p = 0.5; % Momentum (+1,0,-1)
x0_1 = 1.5;
y0_1 = 1.5;
x0_2 = -1.5;
y0_2 = -1.5;
N1 = 1;
N2 = 1;
a = 2;
beta = pi/2;

% Calculated variables
r0 = r0*ones(1,numOfSol);
angles = linspace(0,2*pi*(1-1/numOfSol),numOfSol);
[x0,y0] = pol2cart(angles,r0);
elliptic_x0 = sX*cos(angles);
elliptic_y0 = sY*sin(angles);

% Space definitions
numOfPoints = 256;
samplingSize = 1024; % Might end up being slightly different (bigger)

limitX = 10;
x = linspace(-limitX,limitX,numOfPoints);
dx = x(2)-x(1);
y = x;
dy = dx;
[xx,yy] = meshgrid(x-x0_1,y-y0_1);
r1 = sqrt(xx.^2+yy.^2);
[xx,yy] = meshgrid(x-x0_2,y-y0_2);
r2 = sqrt(xx.^2+yy.^2);
[xx,yy] = meshgrid(x,y);
r = sqrt(xx.^2+yy.^2);

limitZ = 10;
dz = dx^2/8;
sampLeng = floor(length(0:dz:limitZ) / samplingSize);
z = 0:dz*sampLeng:(limitZ+dz*sampLeng);

% Psi
psi = complex(zeros([length(x),length(y),length(z)]));
% Psi_0
switch beamProfile
    case 0
        psi(:,:,1) = exp(-r.^2);
    case 1
        psi(:,:,1) = N1*sech(r);
    case 2
        psi(:,:,1) = N1*sech(r1)+N2*sech(r2); % r1 = sqrt((x-x01).^2+(y-y01).^2)
    case 3
        psi(:,:,1) = sech(r+a).*exp(1i*alpha*(r+a))+sech(r-a).*exp(-1i*alpha*(r-a));
    case 4
        psi(:,:,1) = sech(r+a) + exp(-(r-a).^2);
    case 5
        for n = 1:numOfSol
            psi(:,:,1) = psi(:,:,1) + sech_soliton(x,y,x0(n),y0(n),alphaVec(n),p);
        end
    case 6
        for n = 1:numOfSol
            psi(:,:,1) = psi(:,:,1) + sech_soliton(x,y,elliptic_x0(n),elliptic_y0(n),alphaVec(n),p);
        end
    case 7
        % SPECIAL CASE
        for n = 1:numOfSol
            if n == 2
                p = 1;
                psi(:,:,1) = psi(:,:,1) + sech_soliton(x,y,elliptic_x0(n),elliptic_y0(n),alphaVec(n),p);
            elseif n == 3
                p = -1;
                psi(:,:,1) = psi(:,:,1) + sech_soliton(x,y,elliptic_x0(n),elliptic_y0(n),alphaVec(n),p);
            end
        end
end

psi(:,:,1) = fU;

Kx = linspace(-numOfPoints/2,numOfPoints/2-1,numOfPoints)'*(2*pi/(2*limitX));
Ky = Kx;
[Kxx,Kyy] = meshgrid(Kx,Ky);

fprintf('Finished defining spaces.\n')

%% Propagation

idz = 1i*dz; % Combined factor
Kxy2 = Kxx.^2 + Kyy.^2; % Sum of squared K's
expKxy2 = exp(-idz*Kxy2*0.5); % exp precaulculated for faster calculations
linStep = true; % To choose next step type
tempPsi = complex(zeros(size(psi,1),size(psi,2),2)); % Temporal variable to propagate
tempPsi(:,:,1) = psi(:,:,1);

steps = length(z);
h = waitbar(1/steps,sprintf('Soliton Propagation: %0.02f%%',1/steps*100));

for n = 2:steps
    % Propagate without saving the values
    for m = 1:sampLeng
        % Split-step control, propagate
        if linStep % Lineal step
            tempPsi(:,:,2) = ifft2(ifftshift(expKxy2.*fftshift(fft2(tempPsi(:,:,1)))));
            linStep = false;
        else % Non-lineal step
            tempPsi(:,:,2) = exp(idz*abs(tempPsi(:,:,1)).^2 ./ (1 + 0.02 * abs(tempPsi(:,:,1)).^2)).*tempPsi(:,:,1);
            linStep = true;
        end
        
        % Shift new value into old one
        tempPsi(:,:,1) = tempPsi(:,:,2);
    end
    
    % Save current value
    psi(:,:,n) = tempPsi(:,:,2);
    
    % Update waitbar progress
    if mod(n,10) == 0
        waitbar(n/steps,h,sprintf('Soliton Propagation: %0.02f%%',n/steps*100))
    end
end

fprintf('Finished propagation.\n')
close(h) % Close waitbar

%% Plot
figure(1)
imagesc(x,y,abs(psi(:,:,1)))
title(sprintf('n = %i',1))
colorbar
set(gca,'Ydir','normal')
for n = 2:size(psi,3)
    pause(0.01)
    imagesc(x,y,abs(psi(:,:,n)))
    title(sprintf('n = %i',n))
    colorbar
    set(gca,'Ydir','normal')
end
fprintf('Finished displaying propagation.\n')

if timeScript == 1
    toc
end
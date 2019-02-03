%% Non-Lineal Schrödinger Equation (NLSE)
%
%   dP   1 d^2P   1 d^2P
%  i-- + - ---- + - ---- + |P|^2 P + V P= 0
%   dz   2 dx^2   2 dy^2
%

%%

% Space definitions
numOfPoints = 256;
samplingSize = 512; % Might end up being slightly different (bigger)

windowSize = 10;
x = linspace(-windowSize,windowSize,numOfPoints);
dx = x(2)-x(1);
y = x;
dy = dx;

[xx,yy] = meshgrid(x,y);
r = sqrt(xx.^2+yy.^2);

limitZ = 6;
dz = dx^2/8;
sampLeng = floor(length(0:dz:limitZ) / samplingSize);
z = 0:dz*sampLeng:(limitZ+dz*sampLeng);

% Mode, initial amplitude and omega_0
n = 0;
m = 1;
A0 = 1;
w0 = 1;
[theta,~] = cart2pol(xx,yy);

% Psi
psi = complex(zeros([length(x),length(y),length(z)]));

rw0 = (r/w0);
psi(:,:,1) = A0 * (rw0).^m .* laguerreL(n,2*rw0.^2).^m .* exp(-rw0.^2 + 1i*(m*theta)); % Incomplete
V = -(abs(psi(:,:,1)).^2) + 4*n + 2*m + 2 - 2*r.^2;

Kx = linspace(-numOfPoints/2,numOfPoints/2-1,numOfPoints)'*(2*pi/(2*windowSize));
Ky = Kx;
[Kxx,Kyy] = meshgrid(Kx,Ky);

fprintf('Finished defining spaces.\n')

%% Propagation

idz = 1i*dz; % Combined factor
Kxy2 = Kxx.^2 + Kyy.^2; % Sum of squared K's %%%%%%%%% WHATS WRONG HERE???
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
            tempPsi(:,:,2) = exp(idz*(abs(tempPsi(:,:,1)).^2+V+20)).*tempPsi(:,:,1);
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
for nn = 2:size(psi,3)
    pause(0.01)
    imagesc(x,y,angle(psi(:,:,nn)))
    title(sprintf('n = %i',nn))
    colorbar
    set(gca,'Ydir','normal')
    fprintf('%f\n',sum(sum(abs(psi).^2)))
end
fprintf('Finished displaying propagation.\n')

% Make multiple
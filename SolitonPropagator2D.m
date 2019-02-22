function [psi] = SolitonPropagator2D(psi,V,windowSize,dz,doubleSteps)
    %SOLITONPROPAGATOR2D DESCRIPTION
    %   ETC doubleSteps the number of steps of size 2*dz, or half of steps
    
    % Extract number of points
    numOfPoints = size(psi,1); % Assuming square input
    
    % Calculate K space
    K = linspace(-numOfPoints/2,numOfPoints/2-1,numOfPoints)'*(2*pi/(2*windowSize));
    [Kxx,Kyy] = meshgrid(K,K);
    Kxy2 = Kxx.^2 + Kyy.^2; % Sum of squared K's
    
    % Pre-calculate repeated values
    idz = 1i*dz; % Combined factor
    expKxy2 = exp(-idz*Kxy2*0.5); % exp precaulculated for faster calculations
    
    % Propagator
    for n = 1:doubleSteps
        % Split-step: nonlinear, then linear, two dz propagation
        psi = ifft2(ifftshift(expKxy2.*fftshift(fft2(exp(idz*(abs(psi).^2+V)).*psi))));
    end
end

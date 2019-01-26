function S = sech_soliton(X,Y,X0,Y0,ALPHA,P)
    %SECH_SOLITON Creates a 2D sech at (X0,Y0) within a linspace (X,Y).
    %   The function can also have a phase defined by 2*pi*ALPHA.
    %
    %   SECH_SOLITON(X,Y,X0,Y0) creates the sech centered at (X0,Y0)
    %   without a phase nor momentum. X0 and Y0 are scalars, and X and Y
    %   are vectors.
    %
    %   SECH_SOLITON(X,Y,X0,Y0,ALPHA) adds a phase of 2*pi*ALPHA to the
    %   function.
    %
    %   SECH_SOLITON(X,Y,X0,Y0,ALPHA,P) adds a clockwise or anti-clockwise
    %   momentum to the function, tangent to the vector going from the
    %   origin to (X0,Y0). A positive P gives an anti-clockwise momentum,
    %   while a negative P gives a clockwise momentum.
    
    %   Dependencies: rz_op
    %   https://github.com/jmontemayorm/MATLAB-Functions
    
    %   Author: Javier Montemayor Mancias
    %   Created on: 2018.11.02
    %   Last updated: 2019.01.26
    %   Version: v1.1
    
    % Creates the space for the function, centered in (X0,Y0) and converts
    % to polar.
    [xx,yy] = meshgrid(X-X0,Y-Y0);
    r = sqrt(xx.^2+yy.^2);
    
    % Sets ALPHA to 0 if the argument is not provided.
    if nargin == 4
        ALPHA = 0;
    end
    
    % Sets P to 0 if the argument is not provided (prevents crashing).
    if nargin < 6
        P = 0;
    end
    
    % Creates the function with the respective phase.
    S = sech(r) * exp(2i*pi*ALPHA);
    
    % Adds the momentum.
    if P ~= 0
        m = rz_op(sign(P)*pi/2) * [X0; Y0; 1];
        [xx,yy] = meshgrid(X,Y);
        S = S .* ( exp(abs(P)*1i*m(1)*xx) .* exp(abs(P)*1i*m(2)*yy) );
    end
end
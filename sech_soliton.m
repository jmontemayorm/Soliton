function s = sech_soliton(x,y,x0,y0,alpha,p)
%SECH_SOLITON Creates a 2D sech at (x0,y0) with 2*pi*alpha phase.    
    [xx,yy] = meshgrid(x-x0,y-y0);
    r = sqrt(xx.^2+yy.^2);
    s = sech(r)*exp(2i*pi*alpha);
    if p == 1
        m = rz_op(pi/2)*[x0;y0;1];
        [xx,yy] = meshgrid(x,y);
        %s = s.*(exp(1i*m(1)*xx).*exp(1i*m(2)*yy));
        s = s.*(exp(-3i*xx));
    elseif p == -1
        m = rz_op(-pi/2)*[x0;y0;1];
        [xx,yy] = meshgrid(x,y);
        %s = s.*(exp(1i*m(1)*xx).*exp(1i*m(2)*yy));
        s = s.*(exp(3i*xx));
    end
end
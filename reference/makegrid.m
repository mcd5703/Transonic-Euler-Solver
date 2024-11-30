 function [x,y] = makegrid(NACA,R2,NJ,NK)
    x = zeros(NJ,NK);
    y = zeros(NJ,NK);
    
    % Get the key parameters of the NACA designation
    d1 = floor(NACA/1000);
    m = d1/100;
    NACA = NACA - 1000*d1;
    d2 = floor(NACA/100);
    p = d2/10;
    NACA = NACA - 100*d2;
    t = NACA/100;


    eps = 4*sqrt(3)*t/9;

    z1 = 1;
    z2 = z1 - 4*(1 + 2*eps)/( (1 + 2*eps)^2 + 2*(1 + 2*eps) + 1);

    % Set up airfoil geometry - do a quasi-Newton iteration to make sure
    % that we have constant spacing in the pseudo-circle plane
    
    % Set up target phi space
    phi = linspace(0,2*pi,NK);
    
    % Initial guess at theta
    theta = phi;
    errorMax = 999;
    
    while (errorMax > 1e-10)
        for k = 1:NK
            x(1,k) = 0.5*(1 + cos(theta(k)));
            y(1,k) = 5*t*(.2969*sqrt(x(1,k)) - 0.1260*x(1,k) - 0.3516*x(1,k)^2 + 0.2843*x(1,k)^3 - 0.1036*x(1,k)^4);
            if(theta(k) > pi)
                y(1,k) = -y(1,k);
            end
            if(x(1,k) < p)
                y(1,k) = y(1,k) + m*(2*p*x(1,k) - x(1,k)^2)/p^2;
            else
                y(1,k) = y(1,k) + m*( (1 - 2*p) + 2*p*x(1,k) - x(1,k)^2 )/(1 - p)^2;
            end
            y(1,1) = 0;
            y(1,NK) = 0;
        end

        % Trailing-edge angle
        thetaU = atan2( y(1,2) - y(1,1),x(1,2) - x(1,1) );
        thetaL = atan2( y(1,NK-1) - y(1,NK),x(1,NK-1) - x(1,NK) );
        if (thetaU < 0)
            thetaU = thetaU + 2*pi;
        end
        if (thetaL < 0)
            thetaL = thetaL + 2*pi;
        end
        tauTE = thetaL - thetaU;
        PTE = pi/(2*pi - tauTE);

        % Map the surface
        zeta = zeros(NJ,NK);
        zeta1 = -0.5*(z2 - z1)*PTE;
        zeta2 = 0.5*(z2 - z1)*PTE;
        for k = 1:NK
            RHS = ( (x(1,k) + i*y(1,k) - z1)/(x(1,k) + i*y(1,k) - z2))^PTE;
            zeta(1,k) = (zeta1 - RHS*zeta2)/(1 - RHS);
            if (theta(k) > pi)
                if (imag(zeta(1,k)) > 0)
                    RHS = exp(-2*i*pi*PTE)*RHS;
                    zeta(1,k) = (zeta1 - RHS*zeta2)/(1 - RHS);
                end
            end
        end

        % Find origin
        zeta0 = 0;
        for k = 1:NK-1
            dtheta = theta(k+1) - theta(k);
            zeta0 = zeta0 + (zeta(1,k+1) + zeta(1,k))*dtheta;
        end
        zeta0 = zeta0/(2*pi);
        
        % Measure the angles phi(theta) and associated error
        for k = 1:NK
            phiA(k) = atan2( imag(zeta(1,k)-zeta0), real(zeta(1,k)-zeta0) );
            if (((theta(k)-pi/2)>0)&&(phiA(k) < 0))
                phiA(k) = phiA(k) + 2*pi;
            end
        end
        
        % Set error vectors, phi(theta) - (phi_targ + phi0)
        phi0 = phiA(1);
        phiA(NK) = phi0 + 2*pi;
        error = phiA - (phi + phi0);
        errorMax = max(abs(error));
        
        % Find sensitivities and update
        dtheta = 0.*theta;
        for k = 2:NK-1
            dphidtheta = (phiA(k+1) - phiA(k-1))/(theta(k+1) - theta(k-1));
            dtheta(k) = -error(k)/dphidtheta;
        end
        theta = theta + dtheta;
    end

    % Build pseudo-circular grid
    for k = 1:NK
        R1 = abs(zeta(1,k) - zeta0);
        for j = 1:NJ
            %rad = R1*exp( log(R2/R1)*(j - 1)/(NJ - 1) );
            rad = R1*R2/(R2 + ( (j-1)/(NJ-1) )*(R1 - R2) );
            zeta(j,k) = zeta0 + rad*cos(phiA(k)) + i*rad*sin(phiA(k));
        end
    end

    % Map it
    for j = 1:NJ
        for k = 1:NK
            LHS = ( (zeta(j,k)-zeta1)/(zeta(j,k)-zeta2) )^(1/PTE);
            z = (z2*LHS - z1)/(LHS - 1.);
            x(j,k) = real(z);
            y(j,k) = imag(z);
        end
    end
end


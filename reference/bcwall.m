function [qB] = bcwall(q,qInfPrim,xx,xy,vol,NJ,NK)
    % Free-stream quantities
    rhoInf = qInfPrim(1);
    uInf = qInfPrim(2);
    vInf = qInfPrim(3);
    TInf = qInfPrim(4);
    gamInf = qInfPrim(5);
    
    h0Inf = TInf*gamInf/(gamInf - 1.) + 0.5*(uInf^2 + vInf^2);
    
    % Takes the solution data and zeros-out the normal velocity
    % It extrapolates entropy and fixes total enthalpy
    for j = 1:NJ
        for k = 1:NK
            rho = q(j,k,1);
            u = q(j,k,2)/rho;
            v = q(j,k,3)/rho;
            T = (gamInf-1.)*(q(j,k,4)/rho - 0.5*(u^2 + v^2));
            h0 = q(j,k,4)/rho + T;
            s = T/rho^(gamInf-1.);
            P = rho*T;
            nx = xx(j,k)/sqrt(xx(j,k)^2 + xy(j,k)^2);
            ny = xy(j,k)/sqrt(xx(j,k)^2 + xy(j,k)^2);
            un = u*nx + v*ny;
            ut = u - un*nx;
            vt = v - un*ny;
            
            % Use extrapolated pressure
            %rho = (P/s)^(1/gamInf);
            %T = P/rho;
            
            % Set temperature based on total enthalpy and entropy
            T = (h0 - 0.5*(ut^2 + vt^2))*(gamInf-1.)/gamInf;
            rho = (T/s)^(1./(gamInf-1.));
         
            qB(j,k,1) = rho;
            qB(j,k,2) = rho*ut;
            qB(j,k,3) = rho*vt;
            qB(j,k,4) = rho*(T/(gamInf-1.) + 0.5*(ut^2 + vt^2));
        end
    end
end


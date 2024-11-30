function [R,dt] = resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL)
    % Build the residual, assuming periodicity in K, slip wall at J = 1,
    % and free-stream at J = NJ
    
    % Using a 3rd-order HLLC scheme with Koren limiter

    
    % Some preliminaries
    f = zeros(NJ,NK,4);
    R = zeros(NJ,NK,4);
    rhoInf = qInfPrim(1);
    uInf = qInfPrim(2);
    vInf = qInfPrim(3);
    TInf = qInfPrim(4);
    gamInf = qInfPrim(5);
    h0Inf = TInf*gamInf/(gamInf-1.) + 0.5*(uInf^2 + vInf^2);
    
    sigj = zeros(NJ,NK);
    sigk = zeros(NJ,NK);
    dt = zeros(NJ,NK);
    
    
    % Build spectral radii
    for j = 2:NJ-1
        for k = 1:NK
            rho = q(j,k,1);
            u = q(j,k,2)/rho;
            v = q(j,k,3)/rho;
            T = (gamInf - 1.)*(q(j,k,4)/rho - 0.5*(u^2 + v^2));
            sigj(j,k) = abs(xx(j,k)*u + xy(j,k)*v) + sqrt( (xx(j,k)^2 + xy(j,k)^2)*1.4*T );
            sigk(j,k) = abs(yx(j,k)*u + yy(j,k)*v) + sqrt( (yx(j,k)^2 + yy(j,k)^2)*1.4*T );
            CFLeff = max(CFL/(1 + sqrt(vol(j,k))),1);
            dt(j,k) = CFLeff*vol(j,k)/(max(sigj(j,k),sigk(j,k)));
        end
    end
    
    % Do the K residuals first - again, this is periodic
    for j = 2:NJ-1
        for k = 1:NK-1
            yxA = 0.5*(yx(j,k) + yx(j,k+1));
            yyA = 0.5*(yy(j,k) + yy(j,k+1));
            volA = 0.5*(vol(j,k) + vol(j,k+1));
            
            % Set up useful indices
            km = k-1;
            if (k == 1)
                km = NK-1;
            end
            kp = k+1;
            kpp = k+2;
            if (kpp > NK)
                kpp = 2;
            end
            
            % Get left state
            for n = 1:4
                r = (q(j,kp,n) - q(j,k,n))/(q(j,k,n) - q(j,km,n));
                if ((q(j,k,n) - q(j,km,n)) == 0)
                    r = 1e8;
                end
                phi = max(0, min(2*r, min( (1+2*r)/3, 2)));
                qL(n) = q(j,k,n) + 0.5*phi*(q(j,k,n) - q(j,km,n));
            end
            
            % Get right state
            for n = 1:4
                r = (q(j,k,n) - q(j,kp,n))/(q(j,kp,n) - q(j,kpp,n));
                if ((q(j,kp,n) - q(j,kpp,n)) == 0)
                    r = 1e8;
                end
                phi = max(0, min(2*r, min( (1+2*r)/3, 2)));
                qR(n) = q(j,kp,n) + 0.5*phi*(q(j,kp,n) - q(j,kpp,n));
            end
            
            fhat = HLLCFlux(qL,qR,gamInf,yxA,yyA,volA);
            
            % Subtract free-stream fluxes
            UU = yxA*uInf + yyA*vInf;
            fhat(1) = fhat(1) - rhoInf*UU;
            fhat(2) = fhat(2) - rhoInf*UU*uInf - yxA*rhoInf*TInf;
            fhat(3) = fhat(3) - rhoInf*UU*vInf - yyA*rhoInf*TInf;
            fhat(4) = fhat(4) - rhoInf*UU*h0Inf;
            
            for n = 1:4
                R(j,k,n) = R(j,k,n) + fhat(n);
                R(j,k+1,n) = R(j,k+1,n) - fhat(n);
            end
        end
        R(j,1,:) = R(j,1,:) + R(j,NK,:);
        R(j,NK,:) = R(j,1,:);
    end
    for k = 1:NK
        for j = 1:NJ-1
            xxA = 0.5*(xx(j,k) + xx(j+1,k));
            xyA = 0.5*(xy(j,k) + xy(j+1,k));
            volA = 0.5*(vol(j,k) + vol(j+1,k));
            
            % Set up useful indices
            jm = j-1;
            if (j == 1)
                jm = 1;
            end
            jp = j+1;
            jpp = j+2;
            if (jpp > NJ)
                jpp = NJ;
            end
            
            % Get left state
            for n = 1:4
                r = (q(jp,k,n) - q(j,k,n))/(q(j,k,n) - q(jm,k,n));
                if ((q(jm,k,n) - q(j,k,n)) == 0)
                    r = 1e8;
                end
                phi = max(0, min(2*r, min( (1+2*r)/3, 2)));
                qL(n) = q(j,k,n) + 0.5*phi*(q(j,k,n) - q(jm,k,n));
            end
            
            % Get right state
            for n = 1:4
                r = (q(j,k,n) - q(jp,k,n))/(q(jp,k,n) - q(jpp,k,n));
                if ((q(jp,k,n) - q(jpp,k,n)) == 0)
                    r = 1e8;
                end
                phi = max(0, min(2*r, min( (1+2*r)/3, 2)));
                qR(n) = q(jp,k,n) + 0.5*phi*(q(jp,k,n) - q(jpp,k,n));
            end
            
            fhat = HLLCFlux(qL,qR,gamInf,xxA,xyA,volA);
            
            % Subtract free-stream fluxes
            UU = xxA*uInf + xyA*vInf;
            fhat(1) = fhat(1) - rhoInf*UU;
            fhat(2) = fhat(2) - rhoInf*UU*uInf - xxA*rhoInf*TInf;
            fhat(3) = fhat(3) - rhoInf*UU*vInf - xyA*rhoInf*TInf;
            fhat(4) = fhat(4) - rhoInf*UU*h0Inf;
            
            for n = 1:4
                R(j,k,n) = R(j,k,n) + fhat(n);
                R(j+1,k,n) = R(j+1,k,n) - fhat(n);
            end
        end
    end
    
    % Volume-scale and multiply by DT
    for k = 1:NK
        for j = 2:NJ-1
            R(j,k,:) = dt(j,k)*R(j,k,:)/vol(j,k);
        end
        for n = 1:4
            R(1,k,n) = 0;
            R(NJ,k,n) = 0;
        end
    end
end


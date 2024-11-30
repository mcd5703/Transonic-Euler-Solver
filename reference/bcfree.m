function [qB] = bcfree(q,qInfPrim,xx,xy,vol,NJ,NK)
    % Free-stream quantities
    rhoInf = qInfPrim(1);
    uInf = qInfPrim(2);
    vInf = qInfPrim(3);
    TInf = qInfPrim(4);
    gamInf = qInfPrim(5);

    % Takes the solution data and zeros-out the normal velocity
    for j = 1:NJ
        for k = 1:NK
            % Far-Field BC
            rho = q(j,k,1);
            u = q(j,k,2)/rho;
            v = q(j,k,3)/rho;
            T = (gamInf-1.)*(q(j,k,4)/rho - 0.5*(u^2 + v^2));
            nx = xx(j,k)/sqrt(xx(j,k)^2 + xy(j,k)^2);
            ny = xy(j,k)/sqrt(xx(j,k)^2 + xy(j,k)^2);
            un = u*nx + v*ny;
            unInf = uInf*nx + vInf*ny;
            c = sqrt(gamInf*T);
            cInf = 1;
            ut = u - un*nx;
            vt = v - un*ny;
            utInf = uInf - unInf*nx;
            vtInf = vInf - unInf*ny;
            
            % If we have supersonic inflow, we set to free stream
            if (unInf + cInf < 0)
                qB(j,k,1) = rhoInf;
                qB(j,k,2) = rhoInf*uInf;
                qB(j,k,3) = rhoInf*vInf;
                qB(j,k,4) = rhoInf*(TInf/(gamInf-1.) + 0.5*(uInf^2 + vInf^2));
            else
                % If we have supersonic outflow, we extrapolate
                if (un > c)
                    qB(j,k,:) = q(j,k,:);
                else
                    % Subsonic
                    R1 = un + 2*c/(gamInf - 1.);
                    R2 = unInf - 2*cInf/(gamInf - 1.);
                    % For subsonic outflow
                    if (un > 0)
                        s = T/rho^(gamInf - 1.);
                    % For subsonic inflow
                    else
                        s = TInf/rhoInf^(gamInf - 1.);
                        ut = utInf;
                        vt = vtInf;
                    end
                    un = 0.5*(R1 + R2);
                    c = 0.25*(gamInf - 1.)*(R1 - R2);
                    T = c^2 / gamInf;
                    rho = (T/s)^(1/(gamInf - 1.));
                    u = un*nx + ut;
                    v = un*ny + vt;
                    qB(j,k,1) = rho;
                    qB(j,k,2) = rho*u;
                    qB(j,k,3) = rho*v;
                    qB(j,k,4) = rho*(T/(gamInf - 1.) + 0.5*(u^2 + v^2));
                end
            end
        end
    end
end


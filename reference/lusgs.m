function [dq] = lusgs(R,q,xx,xy,yx,yy,vol,dt,gamInf,NJ,NK)
    dq = zeros(NJ,NK,4);
    dqs = zeros(NJ,NK,4);
    tmp1 = zeros(4,1);
    tmp2 = zeros(4,1);
    sigj = zeros(NJ,NK);
    sigk = zeros(NJ,NK);
    
    % Build spectral radii first
    for j = 1:NJ
        for k = 1:NK
            rho = q(j,k,1);
            u = q(j,k,2)/rho;
            v = q(j,k,3)/rho;
            T = (gamInf - 1.)*(q(j,k,4)/rho - 0.5*(u^2 + v^2));
            sigj(j,k) = abs(xx(j,k)*u + xy(j,k)*v) + sqrt( (xx(j,k)^2 + xy(j,k)^2)*1.4*T );
            sigk(j,k) = abs(yx(j,k)*u + yy(j,k)*v) + sqrt( (yx(j,k)^2 + yy(j,k)^2)*1.4*T );
            sigj(j,k) = sigj(j,k)/vol(j,k);
            sigk(j,k) = sigk(j,k)/vol(j,k);
        end
    end
    % Forward substitution loop
    for j = 1:NJ
        for k = 1:NK
            jm = j-1;
            km = k-1;
            facx = 2;
            facy = 2;
            if (j > 1)
                facx = 1;
                tmp1 = fluxjac(dqs(jm,k,:),q(jm,k,:),xx(jm,k)/vol(jm,k),xy(jm,k)/vol(jm,k),gamInf);
                for n = 1:4
                    tmp1(n) = 0.5*(tmp1(n) + sigj(jm,k)*dqs(jm,k,n));
                end
            end
            if (k > 1)
                facy = 1;
                tmp2 = fluxjac(dqs(j,km,:),q(j,km,:),yx(j,km)/vol(j,km),yy(j,km)/vol(j,km),gamInf);
                for n = 1:4
                    tmp2(n) = 0.5*(tmp2(n) + sigk(j,km)*dqs(j,km,n));
                end
            end
            % Clear out the lower triangular
            D = 1 + dt(j,k)*(facx*sigj(j,k) + facy*sigk(j,k));
            
            for n = 1:4
                dqs(j,k,n) = (-R(j,k,n)*D*vol(j,k) + dt(j,k)*(tmp1(n) + tmp2(n)))/D;
            end
        end
    end
    
    % Clean up periodicity
    for j = 1:NJ
        dqs(j,1,:) = 0.5*(dqs(j,1,:) + dqs(j,NK,:));
        dqs(j,NK,:) = dqs(j,1,:);
    end
    
    % Backward substitution loop
    for j = NJ:-1:1
        for k = NK:-1:1
            jp = j+1;
            kp = k+1;
            facx = 2;
            facy = 2;
            if (j < NJ)
                facx = 1;
                tmp1 = fluxjac(dq(jp,k,:),q(jp,k,:),xx(jp,k)/vol(jp,k),xy(jp,k)/vol(jp,k),gamInf);
                for n = 1:4
                    tmp1(n) = 0.5*(tmp1(n) - sigj(jp,k)*dq(jp,k,n));
                end
            end
            if (k < NK)
                facy = 1;
                tmp2 = fluxjac(dq(j,kp,:),q(j,kp,:),yx(j,kp)/vol(j,kp),yy(j,kp)/vol(j,kp),gamInf);
                for n = 1:4
                    tmp2(n) = 0.5*(tmp2(n) - sigk(j,kp)*dq(j,kp,n));
                end
            end
            % Clear out the upper triangular
            D = 1 + dt(j,k)*(facx*sigj(j,k) + facy*sigk(j,k));
            
            for n = 1:4
                dq(j,k,n) = (dqs(j,k,n) - dt(j,k)*(tmp1(n) + tmp2(n)))/D;
            end
        end
    end
    
    % Clean up periodicity
    for j = 1:NJ
        dq(j,1,:) = 0.5*(dq(j,1,:) + dq(j,NK,:));
        dq(j,NK,:) = dq(j,1,:);
    end
    
    for j = 1:NJ
        for k = 1:NK
            dq(j,k,:) = dq(j,k,:)/vol(j,k);
        end
    end
    %}
end


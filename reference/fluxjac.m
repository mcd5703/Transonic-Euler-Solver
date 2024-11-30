function [fout] = fluxjac(v1,q,xx,xy,gamInf)
    rho = q(1);
    u = q(2)/rho;
    v = q(3)/rho;
    e0 = q(4)/rho;
    phi2 = 0.5*(gamInf - 1)*(u^2 + v^2);
    UU = xx*u + xy*v;
    
    A(1,:) = [0, xx, xy, 0];
    A(2,:) = [xx*phi2 - UU*u, UU - xx*(gamInf-2)*u, xy*u - (gamInf-1)*xx*v, xx*(gamInf-1)];
    A(3,:) = [xy*phi2 - UU*v, xx*v - xy*(gamInf-1)*u, UU - xy*(gamInf-2)*v, xy*(gamInf-1)];
    A(4,:) = [UU*(2*phi2 - gamInf*e0), xx*(gamInf*e0 - phi2) - (gamInf-1)*u*UU, xy*(gamInf*e0 - phi2) - (gamInf-1)*v*UU, gamInf*UU];
    
    fout = zeros(4,1);
    for n = 1:4
        for k = 1:4
            fout(n,1) = fout(n,1) + A(n,k)*v1(k);
        end
    end
end


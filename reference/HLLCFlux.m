function [F] = HLLCFlux(qL,qR,gamInf,xx,xy,vol)
    % Use Roe-averaged wave speeds
    xs = sqrt(xx^2 + xy^2); % Face area
    nx = xx/xs;
    ny = xy/xs;
    
    % Note: Fluxes have volume included, 
    rhoL = qL(1);
    uL   = qL(2)/rhoL;
    vL   = qL(3)/rhoL;
    TL   = (gamInf-1)*(qL(4)/rhoL - 0.5*(uL^2 + vL^2));
    HL   = gamInf*TL/(gamInf-1) + 0.5*(uL^2 + vL^2);
    udotn = uL*nx + vL*ny;
    UUL = udotn;
    FL(1,1) = rhoL*udotn;
    FL(2,1) = rhoL*udotn*uL + rhoL*TL*nx;
    FL(3,1) = rhoL*udotn*vL + rhoL*TL*ny;
    FL(4,1) = rhoL*udotn*HL;
    
    rhoR = qR(1);
    uR   = qR(2)/rhoR;
    vR   = qR(3)/rhoR;
    TR   = (gamInf-1)*(qR(4)/rhoR - 0.5*(uR^2 + vR^2));
    HR   = gamInf*TR/(gamInf-1) + 0.5*(uR^2 + vR^2);
    udotn = uR*nx + vR*ny;
    UUR = udotn;
    FR(1,1) = rhoR*udotn;
    FR(2,1) = rhoR*udotn*uR + rhoR*TR*nx;
    FR(3,1) = rhoR*udotn*vR + rhoR*TR*ny;
    FR(4,1) = rhoR*udotn*HR;
    
    rho  = sqrt(rhoL*rhoR);
    u    = (sqrt(rhoL)*uL + sqrt(rhoR)*uR)/(sqrt(rhoL) + sqrt(rhoR));
    v    = (sqrt(rhoL)*vL + sqrt(rhoR)*vR)/(sqrt(rhoL) + sqrt(rhoR));
    H    = (sqrt(rhoL)*HL + sqrt(rhoR)*HR)/(sqrt(rhoL) + sqrt(rhoR));
    T    = (gamInf-1)*(H - 0.5*(u^2 + v^2))/gamInf;
    c    = sqrt( (gamInf-1)*(H - 0.5*(u^2 + v^2)) );
    
    % Construct eigenvector matrix
    %{
    UU = u*xx + v*xy;
    xs = sqrt(xx^2 + xy^2);
    
    SL = UU - xs*c;
    SR = UU + xs*c;
    SS = UU;
    %SS = ( rhoR*UUR*(SR - UUR) - rhoL*UUL*(SL - UUL) + rhoL*TL - rhoR*TR)/( rhoR*(SR - UUR) - rhoL*(SL - UUL));
    %SS = SS*xs;
    %PLR = 0.5*(rhoL*TL*xs^2 + rhoR*TR*xs^2 + rhoL*(SL - UUL)*(SS - UUL) + rhoR*(SR - UUR)*(SS - UUR));
    
    DS = [0; xx/xs; xy/xs; SS/xs];
    if (SL > 0) % Supersonic flow from left to right
        F = FL;
    elseif (SR < 0) % Supersonic flow from right to left
        F = FR;
    elseif (SS > 0) % Take left intermediate state
        for n = 1:4
            %F(n,1) = ( SS*(SL*qL(n) - FL(n)) + (SL/xs)*(rhoL*TL*xs^2 + rhoL*(SL - UUL)*(SS - UUL))*DS(n))/(SL - SS);
            F(n,1) = ( SS*(SL*qL(n) - FL(n)) + SL*PLR*DS(n)/xs )/(SL - SS);
        end
    else
        for n = 1:4
            %F(n,1) = ( SS*(SR*qR(n) - FR(n)) + (SR/xs)*(rhoR*TR*xs^2 + rhoR*(SR - UUR)*(SS - UUR))*DS(n))/(SR - SS);
            F(n,1) = ( SS*(SR*qR(n) - FR(n)) + SR*PLR*DS(n)/xs )/(SR - SS);
        end
    end
    %}
    
    UU = u*nx + v*ny;
    
    SL = UU - c;
    SR = UU + c;
    %SS = UU;
    %SL = UUL - sqrt(gamInf*TL);
    %SR = UUR + sqrt(gamInf*TR);
    SS = (rhoR*TR - rhoL*TL + rhoL*UUL*(SL - UUL) - rhoR*UUR*(SR - UUR))/( rhoL*(SL-UUL) - rhoR*(SR-UUR));
    
    PLR = 0.5*(rhoL*TL + rhoR*TR + rhoL*(SL - UUL)*(SS - UUL) + rhoR*(SR - UUR)*(SS - UUR));
    DS = [0; nx; ny; SS];
    
    if (SL > 0) % Supersonic flow from left to right
        F = FL;
    elseif (SR < 0) % Supersonic flow from right to left
        F = FR;
    elseif (SS > 0) % Take left intermediate state
        for n = 1:4
            F(n,1) = ( SS*(SL*qL(n) - FL(n)) + SL*(rhoL*TL + rhoL*(SL - UUL)*(SS - UUL))*DS(n))/(SL - SS);
            %F(n,1) = ( SS*(SL*qL(n) - FL(n)) + SL*PLR*DS(n) )/(SL - SS);
        end
    else
        for n = 1:4
            F(n,1) = ( SS*(SR*qR(n) - FR(n)) + SR*(rhoR*TR + rhoR*(SR - UUR)*(SS - UUR))*DS(n))/(SR - SS);
            %F(n,1) = ( SS*(SR*qR(n) - FR(n)) + SR*PLR*DS(n) )/(SR - SS);
        end
    end
    
    % Multiply by face area
    F = xs*F;
end


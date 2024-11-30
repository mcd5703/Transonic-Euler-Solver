clc
clear all
format long

% User-specified simulation parameters
NACA = 2412;
% Freestream Mach number
Minf = 0.82; % CHANGE ME!!!!


%DONT TOUCH
alpha = 0;
gamInf = 1.4;
itmax = 3000;
CFL = 50;
nDisp = 25 ;

%% *****************EXPERT USER SECTION*********************************************************************************************

% User-specified grid parameters
NJ = 41;
NK = 81;
R2 =  100;


% Plotting parameters
set(groot,'DefaultAxesFontSize',20)
set(groot,'DefaultTextInterpreter','Latex')
set(groot,'DefaultLegendInterpreter','Latex')
set(groot,'DefaultAxesTickLabelInterpreter','Latex')

%% DO NOT EDIT BELOW HERE

% Make grid
[x,y] = makegrid(NACA,R2,NJ,NK);
[xx,xy,yx,yy,vol] = metrics(x,y,NJ,NK);

% Initialize solution to free stream
rhoInf = 1; % Density Normalized to 1
uInf = Minf*cos(alpha*pi/180); % Velocity normalized to Mach No.
vInf = Minf*sin(alpha*pi/180);
TInf = 1/gamInf; % Temperature normalized by gamma
qInfPrim = [rhoInf, uInf, vInf, TInf, gamInf];
q = zeros(NJ,NK,4);
for j = 1:NJ
    for k = 1:NK
        q(j,k,1) = rhoInf; %Density
        q(j,k,2) = rhoInf*uInf; %XMomentum
        q(j,k,3) = rhoInf*vInf; %YMomentum
        q(j,k,4) = rhoInf*(TInf/(gamInf-1.) + 0.5*(uInf^2 + vInf^2));
    end
end

% Loop over time
for it = 1:itmax
    %% Calculate update
    
    % Get initial update vector
    q(1,:,:) = 1.5*q(2,:,:) - 0.5*q(3,:,:);
    q(1,:,:) = bcwall(q(1,:,:),qInfPrim,xx(1,:),xy(1,:),vol(1,:),1,NK);
    q(NJ,:,:) = q(NJ-1,:,:);
    q(NJ,:,:) = bcfree(q(NJ,:,:),qInfPrim,xx(NJ,:),xy(NJ,:),vol(NJ,:),1,NK);
    qn = q;
    
    [R,dt] = resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL);
    v1 = lusgs(R,q,xx,xy,yx,yy,vol,dt,gamInf,NJ,NK);
    RN = R;
    
    % Normalize first Krylov vector
    h11 = qdotq(v1,v1,NJ,NK);
    h11 = sqrt(h11);
    v1 = v1/h11;
    
    % Now, get Jacobian-free Arnoldi correction
    %
    % (I + dR/dq)*a1*v1 = -R
    % Minimmize a1*v2 + R
    % a1*v2'*v2 + v2'*R = 0
    %
    %
    q = qn + (2e-8)*v1;
    q(1,:,:) = 1.5*q(2,:,:) - 0.5*q(3,:,:);
    q(1,:,:) = bcwall(q(1,:,:),qInfPrim,xx(1,:),xy(1,:),vol(1,:),1,NK);
    q(NJ,:,:) = q(NJ-1,:,:);
    q(NJ,:,:) = bcfree(q(NJ,:,:),qInfPrim,xx(NJ,:),xy(NJ,:),vol(NJ,:),1,NK);
    [R,dt] = resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL);
    v2 = v1 + (R - RN)/(2e-8);
    
    % Calculate improved update weight
    h22 = qdotq(v2,v2,NJ,NK);
    h12 = qdotq(R,v2,NJ,NK);
    a1 = -h12/h22;
    a1 = max(a1,0.5*h11);
    %}
    % Update, and add BCs back in for good measure
    q = qn + a1*v1;

    if (mod(it,nDisp) == 0) 
        % Plotting data
        rhoresid = norm(norm(R(2:NJ-1,:,1)))/(CFL);
        %fprintf('Iteration %d density residual = %s \n',it,rhoresid)
        for j = 1:NJ
            for k = 1:NK
                Mplot(j,k) = sqrt( (q(j,k,2)^2 + q(j,k,3)^2)/q(j,k,1)^2 / (gamInf*(gamInf-1.)*(q(j,k,4)/q(j,k,1) - 0.5*(q(j,k,2)^2 + q(j,k,3)^2)/q(j,k,1)^2)));
                P(j,k) = (gamInf-1.)*(q(j,k,4) - 0.5*(q(j,k,2)^2 + q(j,k,3)^2)/q(j,k,1));
                Cp(j,k) = (P(j,k) - 1./gamInf)/(0.5*Minf^2); %Normalized using P \gamma
                           end
        end
        
        % Integrate forces on j = 1
        cx = 0;
        cy = 0;
        j = 1;
        for k = 1:NK-1
            cx = cx - 0.5*(Cp(j,k)*xx(j,k) + Cp(j,k+1)*xx(j,k+1));
            cy = cy - 0.5*(Cp(j,k)*xy(j,k) + Cp(j,k+1)*xy(j,k+1));
        end
        cl = cy*cos(alpha*pi/180) - cx*sin(alpha*pi/180);
        cd = cx*cos(alpha*pi/180) + cy*sin(alpha*pi/180);
   
        %fprintf('\t cl = %s, cd = %s \n',cl,cd)
        fprintf('Iter %d: resid = %s, cl = %s, cd = %s \n',it,rhoresid,cl,cd)
%% **********************  END EXPERT USER SECTION  *****************************************************  
       contourf(x,y,Mplot,21)
        xlim([-1 2])
        ylim([-1 1])
        daspect([1 1 1])
        refreshdata
        ylim([-1 1])
        drawnow
    end
end 

Density = q(:,:,1).*0.97;
Pressure = P.*70989.45.*1.4;


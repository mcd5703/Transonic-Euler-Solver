function [xx,xy,yx,yy,vol] = metrics(x,y,NJ,NK)
    % Calculate grid metrics using 2nd-order central differences
    %
    % Hard-coded to assume periodicity in K
    for k = 1:NK
        xxsi(1,k) = x(2,k) - x(1,k);
        yxsi(1,k) = y(2,k) - y(1,k);
        for j = 2:NJ-1
            xxsi(j,k) = 0.5*(x(j+1,k)-x(j-1,k));
            yxsi(j,k) = 0.5*(y(j+1,k)-y(j-1,k));
        end
        xxsi(NJ,k) = x(NJ,k) - x(NJ-1,k);
        yxsi(NJ,k) = y(NJ,k) - y(NJ-1,k);
    end

    for j = 1:NJ
        xeta(j,1) = 0.5*(x(j,2) - x(j,NK-1));
        yeta(j,1) = 0.5*(y(j,2) - y(j,NK-1));
        for k = 2:NK-1
            xeta(j,k) = 0.5*(x(j,k+1)-x(j,k-1));
            yeta(j,k) = 0.5*(y(j,k+1)-y(j,k-1));
        end
        xeta(j,NK) = xeta(j,1);
        yeta(j,NK) = yeta(j,1);
    end

    for j = 1:NJ
        for k = 1:NK
            vol(j,k) = xxsi(j,k)*yeta(j,k) - xeta(j,k)*yxsi(j,k);
            xx(j,k) = yeta(j,k);
            xy(j,k) = -xeta(j,k);
            yx(j,k) = -yxsi(j,k);
            yy(j,k) = xxsi(j,k);
        end
    end
end


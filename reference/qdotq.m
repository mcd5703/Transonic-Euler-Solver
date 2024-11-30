function [H] = qdotq(q1,q2,NJ,NK)
    H = 0;
    for j = 2:NJ-1
        for k = 1:NK
            H = H + q1(j,k,1)*q2(j,k,1);
        end
    end
end


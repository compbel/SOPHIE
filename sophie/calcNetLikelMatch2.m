function likel = calcNetLikelMatch2(AM,degDistrCN,N,enforceTree)

m = sum(sum(AM));
n = length(AM);

if enforceTree == 1
    if m ~= n-1
        likel = -Inf;
        return;
    end
end

d = sum(AM,2)' + sum(AM,1);
flag = 'real';

    C_total = N*degDistrCN;
    C = round(C_total);
    if sum(C) < n
        ind = find(C == 0,n-sum(C),'first');
        C(ind) = 1;
    end
    Delta = find(C >=1,1,'last');
%     N = sum(C(1:Delta));
    D = 1:(N-1);
    M = (D*C_total')/2;

    w = zeros(n,Delta);
    for i = 1:n
        for j = 1:Delta
            if d(i) <= j
                w(i,j) = d(i)*log(j);
%             w(i,j) = 0;
            else
                w(i,j) = -10^9;
            end
        end
    end

    u = cell(1,Delta);
    nAuxFac = 0;
    for j = 1:Delta
        uj = zeros(1,C(j)+1);
        uj(1) = 0;
        for k = 1:(length(uj)-1)
%             uj(k+1) = log(nchoosek(C(j),k));
            uj(k+1) = lognchoosek(C(j),k,flag);
%             uj(k+1) = 0;
        end
        u{j} = uj;
        nAuxFac = nAuxFac + length(uj)-1;
    end
    
    costMatrix = zeros(n,nAuxFac);
    col2deg = zeros(1,nAuxFac);
    for i = 1:n
        c = 1;
        for j = 1:Delta
            for k = 1:(length(u{j})-1)
                costMatrix(i,c) = w(i,j) + (u{j}(k+1) - u{j}(k));
                col2deg(c) = j;
                c = c+1;
            end
        end
    end
    
    costUnmatched = -10^6;
    [Match,uR,uC] = matchpairs(costMatrix,costUnmatched,'max');
    if isempty(uR)
        likel = sum(costMatrix(sub2ind(size(costMatrix), Match(:,1), Match(:,2)))) - m*log(2*M) - lognchoosek(N,n,flag);
    else
        likel = -Inf;
    end
    

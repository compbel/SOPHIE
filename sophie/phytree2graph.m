function [AM,WM,patients] = phytree2graph(Z,n,type)
% AM - adjacency matrix, WM - weight matrix

names = get(Z,'LeafNames');
nSeq = length(names);
dates = zeros(1,2*nSeq-1);
patients = zeros(1,2*nSeq-1);
for i = 1:nSeq
    C = strsplit(names{i},'|');
    dates(i) = str2num(C{3});
    patients(i) = str2num(C{2})+1;
end

[A,idTree,Dist] =  getmatrix(Z);
AM = full(A);
WM = double(AM);
c = 1;
for j = 1:n
    for i = 1:n
        if AM(i,j) > 0
            WM(i,j) = Dist(c);
            c = c + 1;
        end
    end
end

eps = 0.1*min(WM(WM > 0));
outdeg = sum(AM,2);
for i = 1:n
    if (outdeg(i) == 2) && (sum(WM(i,:),2) == 0)
        child = find(AM(i,:));
        if outdeg(child(1)) + outdeg(child(2)) == 0
            if dates(child(1)) > dates(child(2))
                WM(i,child(1)) = eps;
            else
                WM(i,child(2)) = eps;
            end
        end
    end
end

if strcmp(type,'raxml')
    G = digraph(WM);
    root = find(sum(AM,1) == 0);
    D = distances(G,root);
    WM = WM*max(dates)/max(D(~isinf(D)));
end

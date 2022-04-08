function [AMNetCons,W,C] = getConsensusNet(AMNets,likel,nPat)

goodSamp = find(~isinf(likel))';
% optSamp = find(likel == max(likel))';
if isempty(goodSamp)
    AMNetCons = [];
    W = [];
    C = [];
    return;
end
AMNets = AMNets(goodSamp);
likel = likel(goodSamp)';


scaleFactor = -mean(likel);
likel = likel + scaleFactor;
wMatrices = cellfun(@times, AMNets, num2cell(exp(likel)), 'UniformOutput', false);
W=sum(cat(3,wMatrices{:}),3);
C=sum(cat(3,AMNets{:}),3);


[u,v,w] = find(W);
E = [u,v,w];
try
    arb=edmonds(1:nPat,E);
    arbEdgeInd=reconstruct_2(arb); 
catch
    AMNetCons = [];
    W = [];
    C = [];
    return;
end
Earb = E(arbEdgeInd,:);
AMNetCons = zeros(nPat,nPat);
for i = 1:size(Earb,1)
    AMNetCons(Earb(i,1),Earb(i,2)) = 1;
end
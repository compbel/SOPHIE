function [sampTraits,sampLikel,transMatr] = felsenstein_samp_scale_par1_test(AM,WM,traits,Q,traitFreq,nSamp,AMTNtrue,debug)

minP = 1.0e-10;


nNodes = size(AM,1);
nTraits = length(Q);
margDistr = zeros(nNodes,nTraits);
root = find(sum(AM,1) == 0);
parent = zeros(1,nNodes);
for i = 1:nNodes
    if i ~= root
        parent(i) = find(AM(:,i));
    end
end
G = digraph(AM);
Gw = digraph(WM);
dfsorder = flip(dfsearch(G,root))';
outdeg = sum(AM,2)';
leafsCV = (outdeg == 0);
leafs = find(leafsCV);
% traitCounts = histcounts(traits,1:(nTraits+1));


for i = 1:length(leafs)
    margDistr(leafs(i),traits(leafs(i))) = 1;
end

logscaling = zeros(1,nNodes);
transMatr = cell(1,nNodes);
clades = zeros(nNodes,nTraits);
for i = dfsorder
    if outdeg(i) == 0
        clades(i,traits(i)) = 1;
        continue;
    end
    child = find(AM(i,:));
    P1 = expm(Q*WM(i,child(1)));
    P2 = expm(Q*WM(i,child(2)));
    transMatr{child(1)} = P1;
    transMatr{child(2)} = P2;
    clades(i,:) = clades(child(1),:) + clades(child(2),:);
    if (WM(i,child(1)) == 0) && (WM(i,child(2)) == 0)
        margDistr(i,:) = (margDistr(child(1),:) + margDistr(child(2),:))/2;
        continue;
    end
    for j = 1:nTraits
        margDistr(i,j) = (P1(j,:)*(margDistr(child(1),:))')*(P2(j,:)*(margDistr(child(2),:))');
    end
    logscaling(i) = logscaling(child(1)) + logscaling(child(2));
    ind = margDistr(i,:) > 0;
    minLi = min(margDistr(i,ind));
    if minLi < minP
        scFactor = minP/minLi;
        margDistr(i,:) = margDistr(i,:)*scFactor;
        scFactor1 = sum(margDistr(i,:));
        margDistr(i,:) = margDistr(i,:)/scFactor1;
        logscaling(i) = logscaling(i) + log(scFactor) - log(scFactor1);
    end
end


sampTraitsC = cell(1,nNodes);
sampLikel = zeros(nSamp,1);
bfsorder = bfsearch(G,root);
parfor s = 1:nSamp
    s
%     dfsorder_s = flip(dfsorder);
    dfsorder_s = bfsorder';
    sampTrait_s = traits;
    parents_tnet_s = zeros(1,nTraits);
    postProbs = zeros(nNodes,nTraits);
    if debug
        figure
    end
    for i = dfsorder_s
        if leafsCV(i)
%             sampTrait_s(i) = traits(i);
            postProbs(i,:) = margDistr(i,:);
        else
            cladeMult = clades(i,:);
            if i == root
                postProb = traitFreq.*margDistr(i,:);
            else
                par = parent(i);
%                 P = expm(Q*WM(par,i));
                P = transMatr{i};
                postProb = P(sampTrait_s(par),:).*margDistr(i,:);
                noSamp = (parents_tnet_s > 0)&(parents_tnet_s ~= sampTrait_s(par));
                noSamp(sampTrait_s(par)) = 0;
                postProb(noSamp) = 0;
%                 postProb = postProb.*clades(i,:);
                if sum(postProb) == 0
                    sampLikel(s) = -Inf;
                    break;
                end
                cladeMult(sampTrait_s(par)) = max(1,cladeMult(sampTrait_s(par)));
            end
            postProb = postProb.*cladeMult;
            postProb = postProb/sum(postProb,2);
            postProbs(i,:) = postProb;
            CDF = cumsum(postProb);
            p = rand;
            sampTrait_s(i) = find(p <= CDF,1,'first');
%             sampTrait_s(i) = find(postProb == max(postProb),1,'first');
            if (i ~= root) && (sampTrait_s(i) ~= sampTrait_s(par))
                parents_tnet_s(sampTrait_s(i)) = sampTrait_s(par);
            end
        end        
        if i == root
            sampLikel(s) = sampLikel(s) + log(traitFreq(sampTrait_s(i)));
        else
            par = parent(i);
%             P = expm(Q*WM(par,i));
            P = transMatr{i};
            sampLikel(s) = sampLikel(s) + log(P(sampTrait_s(par),sampTrait_s(i)));
        end
        if debug
%             h = plot(Gw,'Layout','layered','NodeLabel',string(sampTrait_s),'EdgeLabel',Gw.Edges.Weight);
            h = plot(G,'Layout','layered','NodeLabel',string(sampTrait_s));
            highlight(h,i,'NodeColor','r','MarkerSize',8);
            [];
        end
    end
    sampTraitsC{s} = sampTrait_s;
    
    if debug
        treeEdgeList = G.Edges;
        treeEdgeList = table2array(treeEdgeList);
        AMNet = getTransNet(treeEdgeList,sampTrait_s,nTraits);
        sens = sum(sum(AMTNtrue.*AMNet,1),2)/sum(sum(AMTNtrue));
%     plotNet(AMNet,AMNet,traits,'traits','');
%     figure
%     plot(G,'Layout','layered','NodeLabel',string(sampTrait_s));
    [];
    end
end

sampTraits = reshape(cell2mat(sampTraitsC),nNodes,nSamp)';
ind = ~isinf(sampLikel);
sampTraits = sampTraits(ind,:);
sampLikel = sampLikel(ind,:);
[sampTraits,ia,ic] = unique(sampTraits,'rows');
sampLikel = sampLikel(ia);

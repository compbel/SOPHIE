function AMNetCons = inferTransNetSOPHIE(filePhylo,nSamp,degDistrType,degDistr,transRate,consType,enforceTree,visualize)

% Input paramers:
% filePhylo - name of file with the phylogenetic tree (in Newick or Nexus formats)
% nSamp - number of samples from the joint distribution of ancestral label assignments
% degDistrType - type of random contact network degree distribution.
%                 possible variants: 'power law', 'custom'
% degDistr - exponent of the power law contact network degree distribution
%         (if degDistrType = 'power law') or actual distribution (if degDistrType = 'custom')
% transRate - transmission rate (used in Jukes-Cantor trait substitution model)
% consType - the likelihood used to build a consensus transmission
%                 tree. Variants: 'joint' (default), 'phylogenetic', 'network'
% enforceTree - whether to force algorithm to consider only sampled
%                acyclic transmission networks. Variants: 0 (default), 1
% visualize -  whether to output a figure with the inferred transmission
%              network. Default value: 0
% 
% Output:
% AMNetCons - adjacency matrix of the inferred transmission network
% 
% Implementation of Edmonds algorithm for maximal arborescence is taken
% from https://www.mathworks.com/matlabcentral/fileexchange/24899-edmonds-algorithm 
% Copyright (c) 2009, Ashish Choudhary. All rights reserved.

treetype = 'favites';
mu = transRate;

tree = phytreeread(filePhylo);
[AMtree,WMtree,patients] = phytree2graph(tree,get(tree,'NumNodes'),treetype);
[AMtree, WMtree,patients] = reduceTree(AMtree, WMtree,patients);


patientList = sort(unique(patients));
patientList = patientList(2:end);
nPat = length(patientList);

qDistr = (1/nPat)*ones(1,nPat)';
[traits,traitFreq] = pat2traits1(patients,patientList);
Q = (mu*(ones(nPat,nPat) - eye(nPat,nPat))).*repmat(qDistr,1,nPat);
Q(1:nPat+1:end) = -sum(Q,2);

tree = digraph(AMtree);
treeEdgeList = tree.Edges;
treeEdgeList = table2array(treeEdgeList);

[sampTraits,felsLikel,~] = felsenstein_samp_scale_par1_test(AMtree,WMtree,traits,Q,traitFreq,nSamp,[],false);

if isempty(sampTraits)
        disp('No networks has been sampled. Check the parameters');
        return;
end

nSampCurr = size(sampTraits,1);
netLikel = zeros(nSampCurr,1);
AMNets = cell(1,nSampCurr);
WMNets = cell(1,nSampCurr);
Ns = zeros(nSampCurr,1);
degDistrsCN = cell(1,nSampCurr);
if strcmp(degDistrType,'power law')
    d = degDistr;
    rzeta = zeta(d);
end

parfor s=1:nSampCurr
    s
    AMNet_s = getTransNet(treeEdgeList,sampTraits(s,:),nPat);
    AMNets{s} = AMNet_s;
    if strcmp(degDistrType,'power law')
        degs = sum(AMNet_s + AMNet_s' > 0);
        degCounts = histcounts(degs,1:nPat);
        Ns(s) = max(ceil(rzeta*degCounts.*((1:(nPat-1)).^d)));
        degDistrsCN{s} = ((1:(Ns(s)-1)).^(-d))/rzeta;
    end
    if strcmp(degDistrType,'custom')
        degDistrsCN{s} = degDistr;
    end
end

parfor s=1:nSampCurr
    s
    AMNet = AMNets{s};
    netLikel(s) = calcNetLikelMatch2(AMNet,degDistrsCN{s},Ns(s),enforceTree);
end

totalLikel = felsLikel + netLikel;
if strcmp(consType,'joint')
    [AMNetCons,~,~] = getConsensusNet(AMNets,totalLikel,nPat);
end
if strcmp(consType,'network')
    [AMNetCons,~,~] = getConsensusNet(AMNets,netLikel,nPat);
end
if strcmp(consType,'phylogenetic')
    [AMNetCons,~,~] = getConsensusNet(AMNets,felsLikel,nPat);
end

if visualize == 1
    plotNet(AMNetCons,[],patientList,'patients','Transmission network');
end


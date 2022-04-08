clear;

filePhylo = 'tree_test.tre';
nSamp = 100000;
degDistrType = 'power law';
degDistr = 2;
transRate = 0.05;
consType = 'joint';
enforceTree = 0;
visualize = 1;

AMNetCons = inferTransNetSOPHIE(filePhylo,nSamp,degDistrType,degDistr,transRate,consType,enforceTree,visualize);
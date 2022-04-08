function [AMtree, WMtree,traits] = reduceTree(AMtree, WMtree,traits)

root = find(sum(AMtree,1) == 0);
tree = digraph(AMtree);
dfsorder = flip(dfsearch(tree,root))';
toKeep = true(1,length(AMtree));
for i = dfsorder
    if (outdegree(tree,i) == 0)
        continue;
    end
    child = successors(tree,i);
    if (traits(child(1)) == traits(child(2))) && (traits(child(1)) > 0) && (traits(child(2)) > 0)
        traits(i) = traits(child(1));
        toKeep(child) = 0;
    end
end
AMtree = AMtree(toKeep,toKeep);
WMtree = WMtree(toKeep,toKeep);
traits = traits(toKeep);
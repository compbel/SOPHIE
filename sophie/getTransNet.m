function AMNet = getTransNet(edgeList,traits,nTraits)
% edgeList = Tree.Edges;
AMNet = zeros(nTraits,nTraits);
for i = 1:size(edgeList,1)
    s = traits(edgeList(i,1));
    t = traits(edgeList(i,2));
    if s ~= t
        AMNet(s,t) = AMNet(s,t) + 1;
    end
end

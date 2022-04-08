function [] = plotNet(AMNetOpt,AMTNtrue,patientList,nodeLabels,figTitle)
        
        nPat = size(AMNetOpt,1);
        GOpt = digraph(AMNetOpt);
        figure
        if strcmp(nodeLabels,'patients')
            h = plot(GOpt,'Layout','force','NodeLabel',patientList);
        end
        if strcmp(nodeLabels,'traits')
            h = plot(GOpt,'Layout','force','NodeLabel',string(1:nPat));
        end
        title(figTitle);
        if ~isempty(AMTNtrue)
            [u,v] = find(AMNetOpt);
            for i = 1:length(u)
                if AMNetOpt(u(i),v(i)) == AMTNtrue(u(i),v(i))
                    highlight(h,u(i),v(i),'EdgeColor','r','LineWidth',1.5);
                end
            end
        end
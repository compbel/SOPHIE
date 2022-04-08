function [] = plotTree(tree,patientList,nodeLabels,traits)
        if strcmp(nodeLabels,'patients')
            figure
            plot(tree,'Layout','layered','NodeLabel',traits2pat(traits,patientList))           
        end
        if strcmp(nodeLabels,'traits')                  
            figure
            plot(tree,'Layout','layered','NodeLabel',string(traits))
        end        
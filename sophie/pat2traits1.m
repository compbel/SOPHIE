function [traits,traitFreq] = pat2traits1(patients,patientList)
traits = zeros(1,length(patients));
for i = 1:length(traits)
    if patients(i) ~= 0
        traits(i) = find(patientList == patients(i));
    end
end
nPat = length(patientList);
traitFreq = (1/nPat)*ones(1,nPat);

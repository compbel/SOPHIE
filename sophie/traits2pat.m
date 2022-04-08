function patients =  traits2pat(traits,patientList)
patients = zeros(1,length(traits));
for i = 1:length(traits)
    patients(i) = patientList(traits(i));
end
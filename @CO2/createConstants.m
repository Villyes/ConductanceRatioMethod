function createConstants()
    clear('functions'); %#ok<CLFUNC>
    file='constants.xls';

    
    a=xlsread(file,'Viscosity','B2:B6');
    save('ViscosityConstants.mat','a');
    
end





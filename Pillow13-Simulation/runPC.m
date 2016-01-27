clear all;
close all;
script0_simulateDataForTesting; % pillow code simulated spikes
        
ind = 1:length(Xsp);
notover = setdiff(ind, tspTot(shortIsis));
Xspover = Xsp;
Xspover(notover,:) = 0;
for i = 1:size(Xspover,2)
    Xspover(notover,i) = 0;
end

noiseless_electrode_data = y0; % noisy_electrode_data = Y;
number_PCs = 3;
Xspnotover = Xsp;
for i = 1:size(Xspnotover,2)
    Xspnotover(tspTot(shortIsis),i) = 0;
end
%pc1 = plot_PC(Xsp,Xspover, Xspnotover, true, false, noiseless_electrode_data, number_PCs);
%pc2 = plot_PC(Xsp,Xspover, Xspnotover,false, true, noiseless_electrode_data, number_PCs);
%pc3 = plot_PC(Xsp,Xspover, Xspnotover, false, false, noiseless_electrode_data, number_PCs);

pc4 = plot_PC1(Xsp,noiseless_electrode_data,number_PCs)
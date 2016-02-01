function pccoord = plot_PC1(X,Y,number_pcs)
%               
%           pccoord = plot_PC(X,Y,number_pcs)
% 
% Calculates the principal components of spike events in raw electrode
% data and finds the projections of each spike into PC space. Plots
% projections into PC space colored by ground truth of which cell causes
% each spike.
% 
%  Input:  
%  ------
%   X [nsamps x ncells] - each column holds spike train of single neuron
%   Y [nsamps x nelec] - raw electrode data
%   number_pcs [1 x 1] - number of principal components to visualize
%
%  Output:  
%  -------
%   pccoord - {nCells}[nspks x number_pcs]
%   plot    - projections of noisy waveforms of each cell onto PCs
%
% 01/20/2016 JRG
%%

% Set number of PCs to calculate
number_pcs_calc = number_pcs;

% Get number of cells and electrodes
nCells = size(X,2);
nElec = size(Y,2);

% Use X vector to find Xall, which records any spike event on any electrode
Xall = find(max(X,[],2)==1);

% For each spike event, capture raw electrode data for 15 bins before and
% 14 bins after (total 30 bins) on all six electrodes.
spikeLen = 60;

gimin = min(find(Xall(1:500)>spikeLen-1));

gimax = max(find(Xall+spikeLen<length(X)));
for gi = gimin:gimax%1:length(Xall)-1;
    sp1temp = [];
    for elecind = 1:nElec
        sp1temp = [sp1temp; Y(Xall(gi)-(spikeLen/2):Xall(gi)+(spikeLen/2-1),elecind)];
    end
    sp1(:,gi) = sp1temp;
end;

% Do PCA on set of all spikes
spCov = sp1*sp1';
% spCov = spCov - mean(spCov(:));
[evec,eval] = eig(spCov);
evals = diag(eval);
[evsort, evind] = sort(evals,'descend');
evecsort = evec(:,evind);

% Find projection of each spike event onto PCs for one cell at a time
% This incorporates "ground truth" of which spikes come from which cell
% and allows projections of each cells' spikes to be plotted in a different
% color.
Xall = sum(X,2); % make spikes single vector to find overlapping
for cellind = 1:nCells
    clear Xg sp2
    Xg = find(X(:,cellind)==1);
    
    % find overlapping spikes
    overlapBins = 30; % max possible is 30
%     minOverlapBins = 25; % min possible is 0
    gimin = min(find(Xg>overlapBins-1));
    gimax = max(find(Xg+overlapBins<length(X)));
    for gi = gimin:gimax
        Xsum(gi) = sum(Xall(Xg(gi)-overlapBins+1:Xg(gi)+overlapBins-1));
%         Xsum(gi) = sum(Xall([Xg(gi)-overlapBins+1:Xg(gi)-minOverlapBins-overlapBins+1 Xg(gi)+minOverlapBins-overlapBins+1:Xg(gi)+overlapBins-1]));

    end
    

    Xg2 = find(Xsum>1);
    % Uncomment in order to plot both overlapping and
    % non-overlapping spikes
%     Xg2 = find(X(:,cellind)==1);
    for gi2 = 1:length(Xg2)-1;        
        sp2temp = [];
        for elecind = 1:nElec
            sp2temp = [sp2temp; Y(Xg(gi2)-(spikeLen/2):Xg(gi2)+(spikeLen/2-1),elecind)];
        end
        sp2(:,gi2) = sp2temp;
    end
    
    for pcind = 1:number_pcs_calc
        
        pccoord{cellind}(:,pcind) = sp2'*evecsort(:,pcind);
    end
    spmedian{cellind} = median(sp2'); % this is the median waveform for each  neuron, used later
    
end


 %% Plot the PCs for each cell together to show clusters
 figure; hold on;
 
 cind = 'brgkmc';
 
 for p1 = 1:3%number_pcs
     for p2 = p1+1:3%number_pcs
         
         for cellind = 1:nCells
             hold on;
             subplot(number_pcs-1,number_pcs-1,(number_pcs-1)*(p1-1)+p2-1);
             hold on; scatter(pccoord{cellind}(:,p1 ),pccoord{cellind}(:,p2),cind(cellind),'x');
             xlabel(sprintf('PC %d',p1)); ylabel(sprintf('PC %d',p2));
             set(gca,'fontsize',14);
         end
         axis equal;
     end
 end
 
 %% Plot the first 3 PCs in 3D
 figure;
 p1 = 1; p2 = 2; p3 = 3;
 for cellind = 1:nCells
     hold on; 
     scatter3(pccoord{cellind}(:,p1 ),pccoord{cellind}(:,p2),pccoord{cellind}(:,p3),cind(cellind),'x');
     xlabel(sprintf('PC %d',p1)); ylabel(sprintf('PC %d',p2)); zlabel(sprintf('PC %d',p3));
     set(gca,'fontsize',14);
     grid on; view(45,30); axis equal;
 end
ph=1;
c1 = median(pccoord{1}); c2 = median(pccoord{2});
hold on; scatter3(c1(1),c1(2),c1(3),50,'filled');

hold on; scatter3(c2(1),c2(2),c2(3),50,'filled')

%%
if elecind == 2
%% Do relative time shifting of waveforms and plot in new figure
% Sort of like a parameterization 

noiseMag = 0.05;
n1 = noiseMag*randn(60,1)';
for cellind = 1:nCells
    spmedianN{cellind} = spmedian{cellind} + noiseMag*randn(120,1)';
end
% basevals = [1 2]; % cell number for base waveform
% shiftvals = [2 1]; % cell number for shifted waveform

basevals =  [1 2 2 3 1 3]; % cell number for base waveform
shiftvals = [2 1 3 1 3 2]; % cell number for shifted waveform

% figure; 
hold on;

for cellind = 1:nCells%*nElec
    clear sumElec sumElecAll
    for elecind = 1:nElec
        % spmedian is of length 120, due to concatenation of 60 bins of the
        % spike on each electrode. Set the middle of each "sum" to the base
        % waveform from each electrode
        sumElec{1} = zeros(180,1); sumElec{1}(61:120,:) = spmedianN{basevals(cellind)}(1:60);
        sumElec{2} = zeros(180,1); sumElec{2}(61:120,:) = spmedianN{basevals(cellind)}(61:120);% + noiseMag*randn(60,1)';
        
        % Shift the other cell's waveform and add it to base waveform
        for tshift = 1:120
            
            sumElecAll{elecind,tshift} = sumElec{elecind};
            sumElecAll{elecind,tshift}(tshift+[1:60]) = ...
                sumElecAll{elecind,tshift}(tshift+[1:60]) + ...
                spmedianN{shiftvals(cellind)}((elecind-1)*60+[1:60])';% + ...; % for second electrode, do time shift
                % noiseMag*randn(60,1);
        end
    end%elecind
    
    % Now that we have all possible shifts of two waveforms on both
    % electrodes, find PCs
    for tshift = 1:120
        for pcind = 1:number_pcs_calc
            % Break out middle 60 bins and of summed waveforms on each
            % electrode and concatenate then project onto PC
            pcshift(tshift,pcind) = ...
                [sumElecAll{1,tshift}(61:120); sumElecAll{2,tshift}(61:120) ]'*evecsort(:,pcind);
        end
        
    end
    % plot PCs

    plot3(pcshift(:,1),pcshift(:,2),pcshift(:,3),cind(cellind),'linewidth',3);
    
    % scatter PCs
    %     for tshift = 1:120
    %         scatter3(pcshift(:,1),pcshift(:,2),pcshift(:,3),cind(cellind));
    %     end

end%cellind
xlabel(sprintf('PC %d',p1)); ylabel(sprintf('PC %d',p2)); zlabel(sprintf('PC %d',p3));
set(gca,'fontsize',14);
grid on; view(45,30); axis equal;
%%
end%if elecind
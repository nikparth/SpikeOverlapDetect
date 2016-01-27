function pccoord = plot_PC(Xsp,Xspover, Xspnotover,over, notover, Y,number_pcs)
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
nCells = size(Xsp,2);
nElec = size(Y,2);

% Use X vector to find Xall, which records any spike event on any electrode


Xall = find(max(Xsp,[],2)==1);

%find the deltas between two cells and for each cell in the overlap, assign
%the same delta value
if over == true
    overcells = zeros(length(Xsp),2);
    deltas = zeros(length(Xsp),1);
    Xall1 = find(max(Xspover,[],2)==1);
 
    for i = 1:2:length(Xall1)
        spike1 = Xall1(i);
        spike2 = Xall1(i+1);
        delta = spike2 - spike1;
       
        
        deltas(spike1) = delta;
        deltas(spike2) = delta;
       
       
        celltmp = [find(Xspover(spike1,:) == 1),find(Xspover(spike2, :) == 1)];
       
        %avg = round(mean(Xall1(i*2-1),Xall1(i*2)));
        overcells(Xall1(i),:) = celltmp;
        overcells(Xall1(i+1),:) = celltmp;
        
        %overcells(avg,:) = celltmp;
        % Xnew(i) = avg;
    end
    % Xall = Xnew;
end


        
    
    
% For each spike event, capture raw electrode data for 15 bins before and
% 14 bins after (total 30 bins) on all six electrodes.
for gi = 1:length(Xall)-1;
    sp1temp = [];
    for elecind = 1:nElec
        sp1temp = [sp1temp; Y(Xall(gi)-15:Xall(gi)+14,elecind)];
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

%if overlapping spikes are being plotted then we use the cellindices
%corresponding to overlaps of two cells (cellind 1 = (2,1) cellind 2 =
%(1,2), etc. (1,2) corresponds to 2 after 1.)
if over == true
    cellInds = zeros(nchoosek(nCells,2)*2,2);
    comb = combntns(1:nCells,2);
    for i = 1:2:size(cellInds)
        cellInds(i:i+1,:) = perms(comb(ceil(i/2.),:));
    end
    for cellind = 1:size(cellInds)
        clear Xg sp2;
        Xg = find(ismember(overcells,cellInds(cellind,:),'rows') == 1);
        length(Xg)
        for gi = 1:length(Xg)-1;        
            sp2temp = [];
            for elecind = 1:nElec
                sp2temp = [sp2temp; Y(Xg(gi)-15:Xg(gi)+14,elecind)];
            end
            sp2(:,gi) = sp2temp;
        end
        length(sp2)
        for pcind = 1:number_pcs_calc
            pccoord{cellind}(:,pcind) = sp2'*evecsort(:,pcind);
        end
       
        
        deltacell{cellind} = deltas(Xg(1:end-1));
    end
     length(pccoord{1})
     length(pccoord{2})
  
else
    for cellind = 1:nCells
        clear Xg sp2
        if over == true
            X = Xspover;
        elseif notover == true
            X = Xspnotover;
        else
            X = Xsp;
        end
        Xg = find(X(:,cellind)==1);

        for gi = 1:length(Xg)-1;        
            sp2temp = [];
            for elecind = 1:nElec
                sp2temp = [sp2temp; Y(Xg(gi)-15:Xg(gi)+14,elecind)];
            end
            sp2(:,gi) = sp2temp;
        end

        for pcind = 1:number_pcs_calc
            
            pccoord{cellind}(:,pcind) = sp2'*evecsort(:,pcind);
          
        end
    end
end

    

 %% Plot the PCs for each cell together to show clusters
 colors = [{'blue'}, {'red'},{'green'}, {'black'}, {'magenta'}, {'yellow'}];
 cind = 'brgkmy';
 figure; hold on;

 for p1 = 1:number_pcs
     
     for p2 = p1+1:number_pcs
         
         T = '\fontsize{16} {';
         for cellind = 1:length(pccoord)
             hold on;
             subplot(number_pcs-1,number_pcs-1,(number_pcs-1)*(p1-1)+p2-1);
             hold on; 
             scatter(pccoord{cellind}(:,p1 ),pccoord{cellind}(:,p2),cind(cellind),'x');
             
             xlabel(sprintf('PC %d',p1)); ylabel(sprintf('PC %d',p2));
             set(gca,'fontsize',14);
             if over == true
                T = [T,'\color{',colors{cellind},'}',num2str(cellInds(cellind,1)),',',num2str(cellInds(cellind,2)), ' '];
                
             else
                 T = [T, '\color{',colors{cellind},'}',num2str(cellind),' '];
             end
         end
         T = strcat(T,'}');
         title(T);
         axis equal;
     end
 end

 
 %% Plot the first 3 PCs in 3D
%  figure;
%  p1 = 1; p2 = 2; p3 = 3;
%  for cellind = 1:length(pccoord)
%      hold on; 
%      scatter3(pccoord{cellind}(:,p1 ),pccoord{cellind}(:,p2),pccoord{cellind}(:,p3),cind(cellind),'x');
%      xlabel(sprintf('PC %d',p1)); ylabel(sprintf('PC %d',p2)); zlabel(sprintf('PC %d',p3));
%      set(gca,'fontsize',14);
%      grid on; view(45,30);
%      if over == true
%              title('\fontsize{16} {\color{blue}1,2 \color{red} 1,3 \color{green} 2,3}');
%              %L = legend('1,2','1,3','2,3', 'location','northeastoutside');
%      else
%          title('\fontsize{16} {\color{blue}1 \color{red} 2 \color{green} 3}');
%      end
%  end

%Plot overlapping spikes in PC space vs. delta overlap
if over == true
    figure;
     p1 = 1; p2 = 2; p3 = 3;
     T = '\fontsize{16} {';
     for cellind = 1:length(pccoord)
         hold on; 
         scatter3(pccoord{cellind}(:,p1 ),pccoord{cellind}(:,p2),deltacell{cellind},cind(cellind),'x');
         xlabel(sprintf('PC %d',p1)); ylabel(sprintf('PC %d',p2)); zlabel(sprintf('spike delta'));
         set(gca,'fontsize',14);
         grid on; view(45,30);
       
         T = [T,'\color{',colors{cellind},'}',num2str(cellInds(cellind,1)),',',num2str(cellInds(cellind,2)), ' '];
         title(strcat(T,' }'));
     end
end
function [fullLifeTotal,fullLifeTotalSE,diffProbTotal,diffRate,diffRateSE] ...
    = plotTransientStatsWrapper(trackDiffDir,fTime,pSize,plotType,nBS)
% PLOTTRANSIENTSTATSWRAPPER is a wrapper function for plotTransientClusterStats
%
% Synopsis: [fullLifeTotal,fullLifeTotalSE,diffProbTotal,diffRate,diffRateSE] ...
%    = plotTransientStatsWrapper(trackDiffDir,fTime,pSize,plotType)
%
%Input
%       trackDiffDir - folder where transient diffusion results are saved
%
%       fTime: time of frame rate
%
%       pSize: pixel size
%
%       plotType: 1 to plot results, 0 otherwise
%
%       nBS: Number of bootstrap samples. Enter 0 to not bootstrap.
%       Input argument added by KJ on 2020-08-24.
%
%
% Output
%       fullLifeTotal: Total fraction of time spent in mobility type from
%       combined data sets
%
%       fullLifeTotalSE: Standard error for fullLifeTotal
%
%       diffProbTotal: Matrix containing the probability of switching from
%       one mobility type to another. Rows indicate the first (original) mobility type,
%       columns indicate the second mobility type. Order of rows: Immobile, confined,
%       free, super. Order of columns: immobile, confined, free, super,
%       no switch, unclassified.
%
%       diffRate: Matrix containing the transition rates between present
%       mobility types. Rows indicate the first(original) mobility type,
%       columns indicate the second mobility type. Order of rows: Immobile, confined,
%       free, super. Order of columns: immobile, confined, free, super.
%
%       diffRateSE: standard error for diffRate.
%
% Tony Vega February 2017
%
% Copyright (C) 2024, Jaqaman Lab - UTSouthwestern 
%
% This file is part of MotionAnalysis_Package.
% 
% MotionAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MotionAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MotionAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Use file directory to load all file names
s = 0; %Used for multiple plots, ignore for now
aux = dir([trackDiffDir,filesep,'*.mat']);
for kk = 1:numel(aux)
    trackDiffCell{kk} = [trackDiffDir filesep aux(kk).name];
end
trackDiffCell= sort(trackDiffCell);

%% Initialize  variables

fractDiff = zeros(length(trackDiffCell),4);

%Median lifetime for each diffusion type for each experiment
lifetimeDistM = NaN(length(trackDiffCell),4);
fullLifeT = NaN(length(trackDiffCell),4);
%Median immobile confinement radius for each experiment
ndRadiiM = NaN(length(trackDiffCell),1);
%Median confined confinement radius for each experiment
cdRadiiM = NaN(length(trackDiffCell),1);
%Median confined diffusion coefficient for each experiment
cdDiffCoefM = NaN(length(trackDiffCell),1);
%Median free diffusion coefficient for each experiment
fdDiffCoefM = NaN(length(trackDiffCell),1);
%Median super diffusion coefficient for each experiment
sdDiffCoefM = NaN(length(trackDiffCell),1);
% Pool all diffusion information
totalTrackDiff = [];

%% Diffusion Information from each experiment

for k = 1:length(trackDiffCell)
    
    % Extract and Store Results
    fileD = load(trackDiffCell{k});
    %Since we don't know what the diffusion data was saved as, we go this
    %roundabout way of changing to a cell so that the name doesn't matter
    cellD = struct2cell(fileD);
    trackDiff = cellD{2};
    if ~isfield(trackDiff,'segmentClass')
        trackDiff = cellD{1};
    end
    totalTrackDiff = [totalTrackDiff; trackDiff];
    [fractDiff(k,:),lifetimeDist,ndRadii,cdRadii,cdDiffCoef,fdDiffCoef,sdDiffCoef,~,fullLifeT(k,:),~] = plotTransientClusterStats(trackDiff,fTime,pSize,0);
    lifetimeDistM(k,:) = nanmedian(lifetimeDist,1);
    ndRadiiM(k) = nanmedian(ndRadii,1);
    cdRadiiM(k) = nanmedian(cdRadii,1);
    cdDiffCoefM(k) = nanmedian(cdDiffCoef);
    fdDiffCoefM(k) = nanmedian(fdDiffCoef);
    sdDiffCoefM(k) = nanmedian(sdDiffCoef);
    
end

%% Total diffusion Information from pooled experiments

[~,lifeTimeDistFull,ndFull,cdFull,cdFullDC,fdFullDC,~,switchingCharTotal,fullLifeTotal,~] = plotTransientClusterStats(totalTrackDiff,fTime,pSize,0);

% Store all confinement radius information
lengthCR = max([length(ndFull),length(cdFull)]);
crMat = NaN(lengthCR,2);
crMat(1:length(ndFull),1) = ndFull;
crMat(1:length(cdFull),2) = cdFull;

% Store all diffusion coefficient information
lengthDC = max([length(fdFullDC),length(cdFullDC)]);
dcMat = NaN(lengthDC,2);
dcMat(1:length(cdFullDC),1) = cdFullDC;
dcMat(1:length(fdFullDC),2) = fdFullDC;

%% Bootstrap data to get standard deviation

if nBS > 0
    
    %Bootstrap-sample information for nBS iterations to get standard deviation
    fullLifeRand = NaN(nBS,4);
    diffRateRand = NaN(4,4,nBS);
    
    for rt = 1:nBS
        
        %Randomly pull with replacement, diffusion tracks for the same number
        %of total diffusion tracks
        trackDiffRand = datasample(totalTrackDiff, length(totalTrackDiff),'Replace',true);
        
        % Get lifetime, switching rates, and ...
        [~,lifetimeRand,~,~,~,~,~,switchingCharRand,fullLifeRand(rt,:),~] = plotTransientClusterStats(trackDiffRand,fTime,pSize,0);
        
        %KJ 200701: Fixing bug in storage of original probabilities below
        diffProb = zeros(4,4);
        diffProbTotal = zeros(4,6);
        
        for dt = 1:4
            
            % Find all switching instances
            test = switchingCharRand(switchingCharRand(:,1)==(dt-1),:);
            % Set transitions to unclassified as 5, no switch has value 6
            test(isnan(test(:,2)),2) = 5;
            [hValues,~] = histcounts(test(:,2),0:6);
            hValues = hValues./sum(hValues);
            index = 1:4;
            nValues = hValues;
            index(dt) = [];
            %Account for the probability of not switching
            nValues(index) = hValues(index).*sum([hValues(index),hValues(6)])./sum(hValues(index));
            %         diffProb(dt,:) = nValues; %KJ 200701: commenting out original code from Tony
            diffProb(dt,:) = nValues(1:4);
            diffProbTotal(dt,:) = hValues;
            
        end
        
        %KJ 200701: commenting out original code from Tony
        %     diffProbTotal = diffProb; %Just in case we ever need everything
        %     diffProb = diffProb(1:4,1:4);
        
        % Collect rates from bootstrap samples
        diffRateRand(:,:,rt) = diffProb./repmat(nanmean(lifetimeRand(:,1:4),1)',1,4);
        
    end
    
    % Get standard deviation of lifetime and transition rate
    diffRateSE = std(diffRateRand,0,3);
    fullLifeTotalSE = std(fullLifeRand);
    
else
    
    diffRateSE = zeros(4,4);
    fullLifeTotalSE = zeros(1,4);
    
end

%% Probability of switching

%KJ 2007071: Fixing bug in storage of original probabilities below
diffProb = zeros(4,4);
diffProbTotal = zeros(4,6);

for dt = 1:4
    
    test = switchingCharTotal(switchingCharTotal(:,1)==(dt-1),:);
    test(isnan(test(:,2)),2) = 5;
    [hValues,~] = histcounts(test(:,2),0:6);
    hValues = hValues./sum(hValues);
    index = 1:4;
    nValues = hValues;
    index(dt) = [];
    nValues(index) = hValues(index).*sum([hValues(index),hValues(6)])./sum(hValues(index));
    %         diffProb(dt,:) = nValues; %KJ 200701: commenting out original code from Tony
    diffProb(dt,:) = nValues(1:4);
    diffProbTotal(dt,:) = hValues;
    
end

%KJ 200701: commenting out original code from Tony
% diffProbTotal = diffProb; %Just in case we ever need everything
% diffProb = diffProb(1:4,1:4);

%% Rate of switching

diffRate = diffProb./repmat(nanmean(lifeTimeDistFull(:,1:4),1)',1,4);

%% Plot for Lifetime, CR, DC

if plotType == 1
    
    figure;
    
    %Fraction of time in diffusion states
    subplot(4,1,1)
    hold on
    scatter(ones(size(fractDiff,1),1)-s,fullLifeT(:,1),'MarkerEdgeColor',[0.5 0.3 0]);
    scatter(2*ones(size(fractDiff,1),1)-s,fullLifeT(:,2),'MarkerEdgeColor','b');
    scatter(3*ones(size(fractDiff,1),1)-s,fullLifeT(:,3),'MarkerEdgeColor','c');
    scatter(4*ones(size(fractDiff,1),1)-s,fullLifeT(:,4),'MarkerEdgeColor','m');
    if size(fractDiff,1)>1
        avgVal = nanmean(fullLifeT(:,1:4),1);
        stdVal = nanstd(fullLifeT(:,1:4),1);
        errorbar((1:4)-s,avgVal,stdVal,'kx','LineWidth',1);
    end
    ylim([0 1])
    set(gca,'XTick',1:4)
    set(gca,'XTickLabel',{'Immobile','Confined','Free','Super'})
    title('Fraction of time spent in diffusion state')
    xlim([0 5])
    
    %Lifetime
    subplot(4,1,2)
    hold on
    title('Median Lifetime')
    scatter(ones(size(fractDiff,1),1)-s,lifetimeDistM(:,1),'MarkerEdgeColor',[0.5 0.3 0]);
    scatter(2*ones(size(fractDiff,1),1)-s,lifetimeDistM(:,2),'MarkerEdgeColor','b');
    scatter(3*ones(size(fractDiff,1),1)-s,lifetimeDistM(:,3),'MarkerEdgeColor','c');
    scatter(4*ones(size(fractDiff,1),1)-s,lifetimeDistM(:,4),'MarkerEdgeColor','m');
    if size(fractDiff,1)>1
        avgVal = nanmean(lifetimeDistM(:,1:4),1);
        stdVal = nanstd(lifetimeDistM(:,1:4),1);
        errorbar((1:4)-s,avgVal,stdVal,'kx','LineWidth',1);
    end
    xlim([0 5])
    set(gca,'XTick',1:4)
    set(gca,'XTickLabel',{'Immobile','Confined','Free','Super'})
    
    %Confinement characteristics
    subplot(4,1,3)
    hold on
    title('Median Confinement')
    scatter(ones(size(fractDiff,1),1)-s,ndRadiiM,'MarkerEdgeColor',[0.5 0.3 0]);
    scatter(2*ones(size(fractDiff,1),1)-s,cdRadiiM,'MarkerEdgeColor','b');
    if size(fractDiff,1)>1
        avgVal = mean([ndRadiiM, cdRadiiM],1);
        stdVal = std([ndRadiiM, cdRadiiM],1);
        errorbar((1:2)-s,avgVal,stdVal,'kx','LineWidth',1);
    end
    xlim([0 3])
    set(gca,'XTick',1:2)
    set(gca,'XTickLabel',{'Immobile','Confined'})
    
    % Diffusion Coefficient
    subplot(4,1,4)%figure;
    hold on
    title('Median Diffusion Coef')
    scatter(ones(size(fractDiff,1),1)-s,cdDiffCoefM,'MarkerEdgeColor','b');
    scatter(2*ones(size(fractDiff,1),1)-s,fdDiffCoefM,'MarkerEdgeColor','c');
    scatter(3*ones(size(fractDiff,1),1)-s,sdDiffCoefM,'MarkerEdgeColor','m');
    if size(fractDiff,1)>1
        avgVal = mean([cdDiffCoefM, fdDiffCoefM,sdDiffCoefM],1);
        stdVal = std([cdDiffCoefM, fdDiffCoefM,sdDiffCoefM],1);
        errorbar((1:3)-s,avgVal,stdVal,'kx','LineWidth',1);
    end
    xlim([0 4])
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',{'Confined','Free','Super'})
    
    %KJ 200824: Save plot
    savefig(fullfile(trackDiffDir,'diffPropertiesFromDC-MSS'));
    
end


end
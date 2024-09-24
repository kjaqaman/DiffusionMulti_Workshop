function [fractDiff,lifetimeDist,ndRadii,cdRadii,cdDiffCoef,fdDiffCoef,...
    sdDiffCoef,switchData,fullLife,diffRate] = plotTransientClusterStats(...
    trackDiff,fTime,pSize,plotType)
%PLOTCLUSTERSTATS extracts a number of characteristics from transient diffusion data
%
% Synopsis: [fractDiff,lifetimeDist,ndRadii,cdRadii,cdDiffCoef,fdDiffCoef,...
%    sdDiffCoef,switchData,fullLife,diffRate] = plotTransientClusterStats(...
%    trackDiff,fTime,pSize,plotType)
%
% INPUT
%       trackDiff: transDiffAnalysisRes output from
%       basicTransientDiffusionAnalysis function
%
%       fTime: time for frame rate
%
%       pSize: pixel size
%
%       plotType: 1 to plot results, 0 otherwise
%
%   REFERENCE
%
%           .momentScalingSpectrum: (Number of classification
%                     subparts)-by-(20+probDim) matrix, where each row
%                     contains the information for one classification
%                     subpart, and the columns store the following:
%                     (1) Start frame of subpart.
%                     (2) End frame of subpart.
%                     (3) Classification of subpart: 1 = confined, 2 =
%                         free, 3 = directed.
%                     (4) MSS slope resulting in classification.
%                     (5-11) Generalized diffusion coefficients for moment
%                            orders 0-6.
%                     (12-18) Scaling power for moment orders 0-6.
%                     (19) Normal diffusion coefficient (from the MSD).
%                     (20) Confinement radius, if subpart is classified as
%                          confined (NaN otherwise).
%                     (21/22/23) Center of subpart, if subpart is
%                                classified as confined (NaN otherwise).
%
% Output
%       fractDiff: 1x4 array of the fraction of the classified tracks
%       exhibiting a type of mobility (immobile, confine, free,super)
%
%       lifetimeDist: array containing the lifetime of segments for a
%       mobility type. Row length is dependent on mobility type with most segments,
%       column length is 4 (immobile, confined, free, super)
%
%       ndRadii: array containg the confinement radii (or rather
%       localization precision) of immobile segments
%
%       cdRadii: array containing the confinement radii of confined
%       segments
%
%       cdDiffCoef,fdDiffCoef,sdDiffCoef: arrays containing the diffusion
%       coefficients of confined, free, and super diffusion tracks
%
%       switchData: 2xnumber of switches, array containing all switching
%       events. First column is first mobility type and second column is
%       the second mobility. For use in plotTransientStatsWrapper,
%       non-switches are also included by having a switch to an artificial
%       state. This is done to get accurate probabilities of switching
%       switching numbers: 0-immobile, 1-confined, 2-free, 3-super,
%       4-artificial state, NaN- unclassified
%
%       fullLife: fraction of the total time(sum of all track durations)
%       spent in a particular mobility type
%
%       diffRate: Matrix containing the transition rates between present
%       mobility types.Rows indicate the first(original) mobility type,
%       columns indicate the second mobility type. Order of rows: Immobile, confined,
%       free, super. Order of columns: immobile, confined, free, super,
%       unclassified, no switch.
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

%% Collect diffusion and switching information

diffData = [];
switchData = [];
count = 0;
for m = 1:size(trackDiff,1)
    for j = 1:size(trackDiff(m).segmentClass,1)
        diffData = [diffData;trackDiff(m).segmentClass(j).momentScalingSpectrum];
        track = trackDiff(m).segmentClass(j).momentScalingSpectrum;
        lastT = size(track,1);%How many segments are there?
        
        if lastT > 1 %If there are multiple, record motion types
            
            switchTypes(:,1) = track(1:lastT-1,3);
            switchTypes(:,2) = track(2:lastT,3);
            switchData = [switchData;switchTypes];
            clear switchTypes
            count = count + 1;
        else
            switchTypes(:,1) = track(1,3);
            switchTypes(:,2) = 4;
            switchData = [switchData;switchTypes];
            clear switchTypes
            
        end
    end
end

%% Update values to reflect real units

diffData(:,20) = (diffData(:,20).*(pSize));
diffData(:,19) = (diffData(:,19).*(pSize^2))./fTime;

%% Diffusion Population

[noDiff,~,~] = find(diffData(:,3) == 0);
[confinedDiff,~,~] = find(diffData(:,3) == 1);
[freeDiff,~,~] = find(diffData(:,3) == 2);
[supDiff,~,~] = find(diffData(:,3) == 3);

fractNo =length(noDiff)/sum(~isnan(diffData(:,3)));
fractFree =length(freeDiff)/sum(~isnan(diffData(:,3)));
fractSup =length(supDiff)/sum(~isnan(diffData(:,3)));
fractCon =length(confinedDiff)/sum(~isnan(diffData(:,3)));
fractDiff = [fractNo,fractCon,fractFree,fractSup];

%% Distribution of life times

lifeTime = (diffData(:,2)-diffData(:,1))+1;
lifeTime = lifeTime.*fTime;
ndLife = lifeTime(noDiff);
cdLife = lifeTime(confinedDiff);
fdLife = lifeTime(freeDiff);
sdLife = lifeTime(supDiff);
ndLifeT = sum(lifeTime(noDiff))./sum(lifeTime(~isnan(diffData(:,3))));
cdLifeT = sum(lifeTime(confinedDiff))./sum(lifeTime(~isnan(diffData(:,3))));
fdLifeT = sum(lifeTime(freeDiff))./sum(lifeTime(~isnan(diffData(:,3))));
sdLifeT = sum(lifeTime(supDiff))./sum(lifeTime(~isnan(diffData(:,3))));
lifetimeDist = NaN(size(diffData,1),4);
lifetimeDist(1:length(ndLife),1) = ndLife;
lifetimeDist(1:length(cdLife),2) = cdLife;
lifetimeDist(1:length(fdLife),3) = fdLife;
lifetimeDist(1:length(sdLife),4) = sdLife;
%     fullLifeT = [lifeTime, diffData(:,3)];
fullLife = NaN(1,4);
fullLife(1,1) = ndLifeT;
fullLife(1,2) = cdLifeT;
fullLife(1,3) = fdLifeT;
fullLife(1,4) = sdLifeT;
%     fullLifeT(1,5) = sum(lifeTime(~isnan(diffData(:,3))));
%{
    labels = {'Immobile','Confined','Free','Super'};
    figure;
    boxplot(fullLife,'notch','on','labels',labels)
%}

%% Distribution of confinement radii

ndRadii = diffData(noDiff,20);
cdRadii = diffData(confinedDiff,20);

%% Distribution of diffusion coefficients

cdDiffCoef = diffData(confinedDiff,19);
fdDiffCoef = diffData(freeDiff,19);
sdDiffCoef = diffData(supDiff,19);


%% Probability of switching

diffProb = zeros(4,6);
xbin = [-0.5 0.5];
for dt = 1:4
    test = switchData(switchData(:,1)==(dt-1),:);
    test(isnan(test(:,2)),2) = 5;
    [hValues,~] = histcounts(test(:,2),0:6);
    hValues = hValues./sum(hValues);
    index = 1:4;
    nValues = hValues;
    index(dt) = [];
    nValues(index) =hValues(index).*sum([hValues(index),hValues(6)])./sum(hValues(index));
    diffProb(dt,:) = nValues;
    xbin = xbin+1;
end
diffProb = diffProb(:,1:4);

% Rate of switching
diffRate = diffProb./repmat(nanmean(lifetimeDist,1)',1,4);

%% Plot for Lifetime, CR, DC

s=0;
if plotType == 1
    
    %Fraction of total time spent in diffusion state
    figure;
    subplot(4,1,1)
    hold on
    scatter(1-s,fullLife(:,1),'MarkerEdgeColor',[0.5 0.3 0]);
    scatter(2-s,fullLife(:,2),'MarkerEdgeColor','b');
    scatter(3-s,fullLife(:,3),'MarkerEdgeColor','c');
    scatter(4-s,fullLife(:,4),'MarkerEdgeColor','m');
    ylim([0 1])
    xlim([0 5])
    set(gca,'XTick',1:4)
    set(gca,'XTickLabel',{'Immobile','Confined','Free','Super'})
    title('Total fraction of time spent in diffusion state')
    
    %Median Lifetime
    subplot(4,1,2)
    hold on
    title('Median Lifetime')
    scatter(ones(size(fractDiff,1),1)-s,nanmedian(lifetimeDist(:,1)),'MarkerEdgeColor',[0.5 0.3 0]);
    scatter(2*ones(size(fractDiff,1),1)-s,nanmedian(lifetimeDist(:,2)),'MarkerEdgeColor','b');
    scatter(3*ones(size(fractDiff,1),1)-s,nanmedian(lifetimeDist(:,3)),'MarkerEdgeColor','c');
    scatter(4*ones(size(fractDiff,1),1)-s,nanmedian(lifetimeDist(:,4)),'MarkerEdgeColor','m');
    medVal = [nanmedian(lifetimeDist(:,1)) nanmedian(lifetimeDist(:,2)) nanmedian(lifetimeDist(:,3)) nanmedian(lifetimeDist(:,4))];
    stdVal = [nanstd(lifetimeDist(:,1)) nanstd(lifetimeDist(:,2)) nanstd(lifetimeDist(:,3)) nanstd(lifetimeDist(:,4))];
    errorbar((1:4)-s,medVal,stdVal,'kx','LineWidth',1);
    xlim([0 5])
    set(gca,'XTick',1:4)
    set(gca,'XTickLabel',{'Immobile','Confined','Free','Super'})
    
    %Confinement characteristics
    subplot(4,1,3)%figure;
    hold on
    title('Median Confinement')
    scatter(ones(size(fractDiff,1),1)-s,nanmedian(ndRadii),'MarkerEdgeColor',[0.5 0.3 0]);
    scatter(2*ones(size(fractDiff,1),1)-s,nanmedian(cdRadii),'MarkerEdgeColor','b');
    medVal = [nanmedian(ndRadii) nanmedian(cdRadii)];
    stdVal = [nanstd(ndRadii) nanstd(cdRadii)];
    errorbar((1:2)-s,medVal,stdVal,'kx','LineWidth',1);
    xlim([0 3])
    %             ylim([0 3])
    %             set(gca,'XTickLabel',{'','','Immobile','','Confined','',''})
    set(gca,'XTick',1:2)
    set(gca,'XTickLabel',{'Immobile','Confined'})
    
    % Diffusion Coefficient
    subplot(4,1,4)%figure;
    hold on
    title('Median Diffusion Coef')
    scatter(ones(size(fractDiff,1),1)-s,nanmedian(cdDiffCoef),'MarkerEdgeColor','b');
    scatter(2*ones(size(fractDiff,1),1)-s,nanmedian(fdDiffCoef),'MarkerEdgeColor','c');
    scatter(3*ones(size(fractDiff,1),1)-s,nanmedian(sdDiffCoef),'MarkerEdgeColor','m');
    medVal = [nanmedian(cdDiffCoef) nanmedian(fdDiffCoef) nanmedian(sdDiffCoef)];
    stdVal = [nanstd(cdDiffCoef) nanstd(fdDiffCoef) nanstd(sdDiffCoef)];
    errorbar((1:3)-s,medVal,stdVal,'kx','LineWidth',1);
    xlim([0 4])
    %             ylim([0 0.5])
    %             set(gca,'XTickLabel',{'','','Confined','','Free','','Super',''})
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',{'Confined','Free','Super'})
    
end

end
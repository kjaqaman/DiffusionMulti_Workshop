classdef TransientDiffusionAnalysisProcess < PostTrackingProcess
    % A concrete class for analyzing transient track diffusion
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
    
    

    % Tony Vega, January 2018, Adapted from MotionAnalysisProcess
    
    methods (Access = public)
        function obj = TransientDiffusionAnalysisProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                super_args{1} = owner;
                super_args{2} = TransientDiffusionAnalysisProcess.getName;
                super_args{3} = @analyzeTransientMovieMotion;
                if isempty(funParams)  % Default funParams
                    funParams = TransientDiffusionAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@PostTrackingProcess(super_args{:});
        end

        function h=draw(obj,iChan,varargin)
            h = obj.draw@PostTrackingProcess(iChan,varargin{:},'useCache',true);
        end

        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'diffAnalysisRes', 'tracks'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',[],@(x) isempty(x) || isscalar(x) && obj.checkFrameNum(x));
            ip.addParamValue('iZ',[], @(x)ismember(x,1:obj.owner_.zSize_));
            ip.addParamValue('useCache',false,@islogical);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            nOutput = numel(output);
            
            % Data loading

            % load outFilePaths_{1,iChan}
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, output{:});
            
            varargout = cell(nOutput);
            for i = 1:nOutput
                switch output{i}
                    case 'tracks'
                        tracksFinal = s.(output{i});
                        if ~isempty(iFrame),
                            % Filter tracks existing in input frame
                            trackSEL=getTrackSEL(tracksFinal);
                            validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                            [tracksFinal(~validTracks).tracksCoordAmpCG]=deal([]);
                            
                            for j=find(validTracks)'
                                tracksFinal(j).tracksCoordAmpCG = tracksFinal(j).tracksCoordAmpCG(:,1:8*(iFrame-trackSEL(j,1)+1));
                                if isempty( obj.displayMethod_)
                                    tracksFinal(j).framesUsed =1:(iFrame-trackSEL(j,1)+1);
                                else
                                    tracksFinal(j).framesUsed =max([1,(iFrame-trackSEL(j,1)+1)-obj.displayMethod_{1,1}.dragtailLength]):iFrame;
                                end
                            end
                            varargout{i} = tracksFinal;
                        else
                            varargout{i} = tracksFinal;
                        end
                    case 'diffAnalysisRes'
                        varargout{i} = s.(output{i});
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            types = TransientDiffusionAnalysisProcess.getTrackTypes();
            colors = vertcat(types.color);
            
            output(1).name='Classified tracks';
            output(1).var='tracks';
            output(1).formatData=@TransientDiffusionAnalysisProcess.formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)TracksDisplay('Color', colors);
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Transient Motion Analysis';
        end
        function h = GUI()
            h = @transientDiffusionAnalysisProcessGUI;
        end
     
        function alpha = getAlphaValues()
            alpha=[0.01 0.05 0.1 0.2];
        end
        
        %Unclear whether this is necessary...
        function peakA = getPeakAlpha()
            peakA=95;
        end
        
        function methods = getConfinementRadiusMethods()
            methods(1).type = 0;
            methods(1).name = 'Mean positional standard deviation';
            methods(2).type = 1;
            methods(2).name = 'Minimum positional standard deviation';
%             methods(3).type = 2;
%             methods(3).name = 'Rectangle approximation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'TransientDiffusionAnalysis'];

            if owner.is3D
                funParams.probDim = 3;
            else
                funParams.probDim = 2;
            end
            
            %funParams.checkAsym = 0; %Not functional yet
            %funParams.alphaAsym = 0.05;
            funParams.plotRes = 0;
            funParams.peakAlpha = 95;
        end
        
        function displayTracks = formatTracks(tracks)
            % Format classified tracks into structure for display
            
            % Read track types and classification matrix
            types = TransientDiffusionAnalysisProcess.getTrackTypes();
            % Matrix of the number of track segments by 3
            track_class = vertcat(tracks.classification);

            % Number of labels per track needed
            nLabels = cellfun('size',{tracks.classification},1);
            % Labels is a cell array of indices corresponding to to track segment types,
            %   see getTracksTypes
            % Initialize all labels as unlabeled
            labels = arrayfun(@(x) ones(x,1)*numel(types),nLabels,'UniformOutput',false);
            % Map indicates position of last label for each track
            map = cumsum(nLabels);
            
            % logical array per track index of if there is more than one label per track
            nLabels_gt_1 = nLabels > 1;
            % labels index where the labels for each track starts
            %  if there is more than one label per track
            start = map(nLabels_gt_1)-nLabels(nLabels_gt_1)+1;
            % labels index where the labels for each track ends
            %  if there is more than one label per track
            finish = map(nLabels_gt_1);

            % Set labels as needed
            for i = 1 : numel(types) - 1
                % idx is a logical array if each track _segment_ belongs to types(i)
                idx = types(i).f(track_class);
                % idx2 is a cell array of logical arrays
                %  the number of cells corresponds to each track index
                %  the index of the logical array in each cell refers to each track segment
                % Setup logical arrays. This works for when nLabels == 1
                idx2 = num2cell(idx(map));
                % deal with nLabels > 1 separately, grab range of labels for each track segment
                %  corresponding to each track        
                idx2(nLabels_gt_1) = arrayfun(@(s,e) idx(s:e),start,finish,'UniformOutput',false);
                % idx is now a cell array of logical arrays marking for each track
                %  if each track segment belongs to types(i)
                idx = idx2;
                % Select only the indices where at least one segment belongs to type(i)
                any_idx = cellfun(@any,idx);
                % Assign label as i for each track segment belonging to types(i)
                % Merge with previous labels
                labels(any_idx) = cellfun(@(idx_i,labels_i) labels_i.*~idx_i + i.*idx_i, ...
                    idx(any_idx), labels(any_idx)','UniformOutput',false);
            end
            % Assign labels to each track
            [tracks.label] = deal(labels{:});
            %Count for storing extra segments
            count = 1;
            trackAmendIdx = [];
            %Go through tracks with multiple labels
            for m =  find(nLabels_gt_1)
                useClass = zeros(size(tracks(m).transientClassification,1),1);
                for k = 1:size(tracks(m).transientClassification,1)
                    if intersect(tracks(m).transientClassification(k,1):tracks(m).transientClassification(k,2),tracks(m).framesUsed)
                        useClass(k) = 1;
                    end
                end
                %Find out which tracks or track segments are in a given
                %time frame
                if sum(useClass) ==1
                    %If the number of unique tracks is the same as the
                    %number of classifications, and only one is present in
                    %the time frame, no formatting is needed
                    if length(unique(tracks(m).classification(:,4)))-size(tracks(m).classification,1) == 0 && size(tracks(m).tracksFeatIndxCG,1) >1
                        continue
                    % If there are less unique tracks than classification
                    % (i.e transient diffusion), only retain full tracks
                    % (whether present or not) and track segments that are
                    % present
                    elseif length(unique(tracks(m).classification(:,4)))-size(tracks(m).classification,1) < 0 && size(tracks(m).tracksFeatIndxCG,1) >1
                        idx = find(useClass);
                        holdClass = tracks(m).classification(:,4);
                        holdClass(idx) = holdClass(idx)+rand(length(idx),1);%temporarily change this value so it's not found in the next step
                        [~,idxUniq,~]=unique(holdClass,'stable');
                        if ~isempty(idxUniq)
                            holdClass(~ismember(1:length(holdClass),idxUniq)) = tracks(m).classification(idx,4);
                        end
                        idxKeep = ~ismember(holdClass,tracks(m).classification(idx,4));
                         
                        tracks(m).label =tracks(m).label(idxKeep);
                    % If there's only one track present, keep the relevant
                    % segment
                    else
                        tracks(m).label = tracks(m).label(useClass==1);
                    end
                    
                %If track doesn't exist in given frame range, for now
                %just give it a single arbitrary value to avoid code
                %failure    
                elseif sum(useClass) ==0
                    if length(unique(tracks(m).classification(:,4)))-size(tracks(m).classification,1) == 0
                        continue
                    else
                        tracks(m).label = tracks(m).label(1);
                    end
                % More complicated scenarios    
                else
                    %If many non-transient tracks are simply
                    %present, no formatting is needed
                    if length(unique(tracks(m).classification(:,4)))-size(tracks(m).classification,1) == 0 
                        continue
                    % If there are less unique tracks than classification
                    % (i.e transient diffusion), only retain full tracks
                    % (whether present or not) and track segments that are
                    % present. In this case, many tracks/segments are present   
                    elseif length(unique(tracks(m).classification(:,4)))-size(tracks(m).classification,1) < 0 %&& size(tracks(m).tracksFeatIndxCG,1) >1
                        idx = find(useClass);
                        
                        %If it is the simple case of one segment per track
                        %is present in a time frame
                        if length(unique(tracks(m).classification(useClass==1,4))) == length(idx)
                            
                            holdClass = tracks(m).classification(:,4);
                            holdClass(idx) = holdClass(idx)+rand(length(idx),1);%temporarily change this value so it's not found in the next step
                            [~,idxUniq,~]=unique(holdClass,'stable');
                            if ~isempty(idxUniq)
                                holdClass(~ismember(1:length(holdClass),idxUniq)) = tracks(m).classification(idx(1),4);
                            end
                            idxKeep = ~ismember(holdClass,tracks(m).classification(idx,4));

                            tracks(m).label =tracks(m).label(idxKeep);
                        %The most involved scenario, many segments from one
                        %or more tracks are present in a given time frame
                        else
                            %And here we do the painful job of separating track information
                            %Hold onto label information before changing things
                            holdLabels = tracks(m).label;
                            %Get only labels that are present in time frame
                            holdLabelsN = holdLabels(idx);
                            %Perform analysis on a track by track basis
                            %Start replicating track information for other segments, don't
                            %change first segment yet!
                            [~,iS,~]=unique(tracks(m).classification(useClass==1,4),'stable');
                            [~,iSFull,~]=unique(tracks(m).classification(:,4),'stable'); %Verify
                            exTrackIdx = find(~ismember(idx,idx(iS))); %Be sure this doesn't find irrelevant segments
                            tracks(m).label = holdLabels(iSFull);

                            for n = exTrackIdx'
                                trackIdx = tracks(m).transientClassification(idx(n),4);
                                trackSeg = tracks(m);
                                segClass = tracks(m).transientClassification(idx(n),1:2);
                                %For a second, let's pretend all tracks are
                                %present. Defeats the purpose obviously but
                                %may isolate error.
                                trackSeg.tracksCoordAmpCG = tracks(m).tracksCoordAmpCG(trackIdx,segClass(1)*8-7:min([length(tracks(m).tracksCoordAmpCG),segClass(2)*8]));
%                                 trackSeg.tracksCoordAmpCG = tracks(m).tracksCoordAmpCG(trackIdx,:);
                                trackSeg.label = holdLabelsN(n);
                                trackSeg.tracksFeatIndxCG = tracks(m).tracksFeatIndxCG(trackIdx,:);
                                trackSeg.seqOfEvents = tracks(m).seqOfEvents([1,end],:);
                                trackSeg.seqOfEvents(2,3) =1;
                                trackAmend(count) = trackSeg;
                                trackAmendIdx = [trackAmendIdx;m];
                                count = count+1;
                            end

                            %Okay now we change the original track info
                            segClass = tracks(m).transientClassification(idx(iS),1:2);
                            tracksCoordAmpCGMat = NaN(size(tracks(m).tracksFeatIndxCG,1),size(tracks(m).tracksCoordAmpCG,2));
                            tracksFeatIndxCGMat = NaN(size(tracks(m).tracksFeatIndxCG,1),size(tracks(m).tracksFeatIndxCG,2));
                            for k = 1:size(segClass,1)
                                trackIdx =tracks(m).transientClassification(idx(iS(k)),4);
                                %For a second, let's pretend all tracks are
                                %present. Defeats the purpose obviously but
                                %may isolate error.
%                                 tracksCoordAmpCGMat(k,:) = tracks(m).tracksCoordAmpCG(trackIdx,:);
                                tracksCoordAmpCGMat(trackIdx,segClass(k,1)*8-7:min([length(tracks(m).tracksCoordAmpCG),segClass(k,2)*8])) = tracks(m).tracksCoordAmpCG(trackIdx,segClass(k,1)*8-7:min([length(tracks(m).tracksCoordAmpCG),segClass(k,2)*8]));
                                tracksFeatIndxCGMat(trackIdx,:) = tracks(m).tracksFeatIndxCG(trackIdx,:);
                            end
% %                             msIdx = ismember(tracks(m).seqOfEvents(:,3),tracks(m).transientClassification(idx(iS),4));
% %                             seqOfEventsMat = tracks(m).seqOfEvents(msIdx,:);
% %                             [~,~,c] = unique(seqOfEventsMat(:,3),'stable');
% %                             seqOfEventsMat(:,3) = c;
% %                             tracks(m).seqOfEvents =seqOfEventsMat;
                            tracks(m).tracksCoordAmpCG =tracksCoordAmpCGMat;
                            tracks(m).tracksFeatIndxCG =tracksFeatIndxCGMat;
                            
                            %{
                            %Set the original to the one of the classifications
                            tracks(m).label = holdLabels(1);

                            %Start replicating track information for other segments, don't
                            %change first segment yet!
                            for n = 2:size(holdLabels,1)
                                trackSeg = tracks(m);
                                segClass = tracks(m).transientClassification(idx(n),1:2);
                                trackSeg.tracksCoordAmpCG = tracks(m).tracksCoordAmpCG(segClass(1)*8-7:min([length(tracks(m).tracksCoordAmpCG),segClass(2)*8]));
                                trackSeg.label = holdLabels(n);
                                trackAmend(count) = trackSeg;
                                trackAmendIdx = [trackAmendIdx;m];
                                count = count+1;
                            end
                            %Okay now we change the original track info
                            segClass = tracks(m).transientClassification(idx(1),1:2);
                            tracks(m).tracksCoordAmpCG = tracks(m).tracksCoordAmpCG(segClass(1)*8-7:min([length(tracks(m).tracksCoordAmpCG),segClass(2)*8]));
                            %}
                        end
                    end

                end
            end
            
            %If there were additional segments, attach to end
            %Keep old tracks
            if count ==0
                oldTracks = tracks;
                tracks(size(tracks,1)+1:size(tracks,1)+size(trackAmend,2)) = trackAmend;
            end
            
            % Format tracks using TrackingProcess utility function
            displayTracks = TrackingProcess.formatTracks(tracks);
            
            %After formatting tracks, correct track indices for additional
            %segments. Keep in mind this may change with M/S
            if count >1
                displayTracksA = TrackingProcess.formatTracks(trackAmend);
                for n = 1:length(trackAmendIdx)
                    displayTracksA(n).number = trackAmendIdx(n);
                end
               displayTracks(size(displayTracks,1)+1:size(displayTracks,1)+size(trackAmend,2)) = displayTracksA; 
            end
        end
        
        function types = getTrackTypes()
            % Get the color map for classified tracks
            %
            % see also: plotTracksDiffAnalysis2D
            % Immobile: brown
            types(1).name = 'immobile';
            types(1).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 0;
            types(1).color = [0.5 0.3 0];
            % Linear 1D confined: orange
            types(2).name = 'linear & 1D confined diffusion';
            types(2).f = @(x) x(:, 1) == 1 & x(:, 3) == 1;
            types(2).color = [1 0.7 0];
            % Linear 1D normal: bright red 
            types(3).name = 'linear & 1D normal diffusion';
            types(3).f = @(x) x(:, 1) == 1 & x(:, 3) == 2;
            types(3).color = [1 0 0];
            % Linear 1D super: bright green
            types(4).name = 'linear & 1D super diffusion';
            types(4).f = @(x) x(:, 1) == 1 & x(:, 3) == 3;
            types(4).color = [0 1 0];
            % Linear 1D too short: yellow
            types(5).name = 'linear & too short to analyze 1D diffusion';
            types(5).f = @(x) x(:, 1) == 1 & isnan(x(:, 3));
            types(5).color = [1 1 0];
            % Random/unclassified & 2D confined: blue
            types(6).name = 'random/unclassified & 2D confined diffusion';
            types(6).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 1;
            types(6).color = [0 0 1];
            % Random/unclassified & 2D normal: cyan
            types(7).name = 'random/unclassified & 2D normal diffusion';
            types(7).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 2;
            types(7).color = [0 1 1];
            % Random/unclassified & 2D super: magenta
            types(8).name = 'random/unclassified & 2D super diffusion';
            types(8).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 3;
            types(8).color = [1 0 1];
            % Random & 2D too short: purple
            types(9).name = 'random & too short to analyze 2D diffusion';
            types(9).f = @(x) x(:, 1) == 0 & isnan(x(:, 2));
            types(9).color = [.6 0 1];
            % Too short: grey
            types(10).name = 'too short for any analysis';
            types(10).f = @(x) 1;
            types(10).color = [.7 .7 .7];
        end
        
        function hIm = showTrackTypes(newFigure)
            % MotionAnalysisProcess.showTypes: Show types and their colors
            
            if(nargin < 1)
                newFigure = false;
            end
            if(newFigure)
                figure;
            end
            
            types = TransientDiffusionAnalysisProcess.getTrackTypes;
            
            % Create a grey background
            hIm = imshow(ones(300)*0.3);
            % Display type names in their color from top to bottom
            for i=1:length(types)
                text(1,length(types)-i+1,types(i).name, ...
                    'Color',types(i).color, ...
                    'Units','characters');
            end
        end
        
    end
end

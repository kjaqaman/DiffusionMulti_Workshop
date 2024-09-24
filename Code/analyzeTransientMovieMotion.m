function analyzeTransientMovieMotion(movieDataOrProcess,varargin)
% analyzeTransientMovieMotion detects transient changes in particle tracks
%
% SYNOPSIS analyzeTransientMovieMotion(movieDataOrProcess,varargin)
%
% INPUT
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       A character string specifying the directory where to save
%       the motion analysis results.
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       The integer index of the channel(s) containing tracks to be analyzed.
%
%       ('probDim' -> Positive integer)
%       Problem dimensionality. Default is 2.
%
%       ('peakAlpha' -> Positive integer)
%       Confidence level (in percent) for choosing peaks when initially
%       segmenting track. Default: 95.
%
%       ('plotRes' -> Boolean)
%       Whether to plot results. Default:0.
%
% Tony Vega, Jan 2018 , Adapted from analyzeMovieMotion
% Updated by KJ, May 2019
%
% Jesus Vega-Lugo, September 2023, Updated to allow the used of clean
% tracks coming from refineTracksWithinMask
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

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get the MotionAnalysisProcess and create it if it does not exist
[movieData, postProc] = getOwnerAndProcess(movieDataOrProcess,'TransientDiffusionAnalysisProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(postProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check tracking process first
iTrackProc = movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to post-processing!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  post-processing!']);

%load clean tracks
if isfield(p,'useCleanTracks') && p.useCleanTracks
    iTrackProc = movieData.getProcessIndex('RefineTracksWithinMaskProcess',1);

    assert(~isempty(iTrackProc),['No cleaned tracks were found! '...
    'Please run RefineTracksWithinMaskProcess to get clean tracks!'])
    
    trackProc = movieData.getProcess(iTrackProc);
end

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
end
postProc.setInFilePaths(inFilePaths);

% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory);
postProc.setOutFilePaths(outFilePaths);

%% --------------- transient diffusion analysis ---------------%%%

disp('Begin analyzing track diffusion ...')

for i = p.ChannelIndex
    
    tracks = trackProc.loadChannelOutput(i);
    diffAnalysisRes = basicTransientDiffusionAnalysisv1(tracks,...
        p.probDim,p.plotRes,p.peakAlpha);
    
    for j = 1:numel(tracks)
        fullClassHandler = [];
        fullTransInfo = [];
        for k = 1:size(tracks(j).tracksFeatIndxCG,1)
            classHandler = NaN(size(diffAnalysisRes(j).segmentClass(k).momentScalingSpectrum,1),4);
            classHandler(1:size(classHandler,1),2) =diffAnalysisRes(j).segmentClass(k).momentScalingSpectrum(:,3);
            classHandler(1:size(classHandler,1),4) =k;
            transInfo = diffAnalysisRes(j).segmentClass(k).momentScalingSpectrum(:,1:3);
            transInfo(:,4) =k;
            fullClassHandler = [fullClassHandler;classHandler];
            fullTransInfo = [fullTransInfo;transInfo];
        end
        tracks(j).classification = fullClassHandler;
        tracks(j).transientClassification = fullTransInfo;
    end
    
    % save each projData in its own directory
    save(outFilePaths{1,i},'diffAnalysisRes', 'tracks')
    
end

disp('Finished analyzing track diffusion!')

 function analyzeDiffusionModesMLMD(MLMD,varargin)
%ANALYZEDIFFUSIONMODESMLMD performs diffusion mode analysis on single particle tracks
%
%SYNOPSIS analyzeDiffusionModesMLMD(MLMD,varargin)
%
%INPUT 
%   MLMD: MovieList or MovieData object for movie(s) to be analyzed
%
%OPTIONAL (as name-value pairs)
%
%   channel2analyze: integer index of the channel(s) containing tracks to be analyzed
%
%         minLength: Minimum length of a track to be included in analysis. 
%                    Default: 5 frames.
%
%             alpha: Alpha-value for the statistical test to determine number of modes. 
%                    Default: 0.01.
%
%        maxNumMode: Upper limit on the number of modes. 
%                    Default: 10.
%
%       binStrategy: Binning strategy for calculating the cumulative histogram. 
%                    1 for using "histogram" and 2 for using the data directly. 
%                    Default: 2.
%
%       subSampSize: size of subsample to use in mode decomposition. In this 
%                    case, the original data are subsampled many times, each 
%                    time with the specified subSampSize. The output is then 
%                    the diffusion mode decomposition result for each of the 
%                    subsamples. Enter [] if no sub-sampling. 
%                    Default: [].
%
%         doControl: 1 in order to do mono-exponential control, 0 otherwise. 
%                    Default: 1.
%
%    forceDecompose: 1 in order to force run diffusion mode decomposition, 0 
%                    to use old decomposition results if available.  
%                    Default: 1.
%
%          showPlot: 1 to plot the histogram of square displacements and  
%                    fitted exponentials. 0 to not plot anything. 
%                    Default: 1.
%
%          plotName: string with the title of the plotted figure. 
%                    Default: 'Figure'.
%
%    useCleanTracks: 1 to use cleaned tracks, 0 to use original tracks
%                    Default: 0
%
%   PARAMETERS FOR TRACK DIFFUSION MODE ASSIGNMENT
%
%   diffModeDividerStruct: See function trackDiffModeAnalysis for structure details.       
%                          Default: [], in which case the tracks are not
%                          classified into modes but only their diffusion
%                          coefficient is output.
%
%OUPUT 
%
%      Output is saved in directory DiffusionModeAnalysis belonging to each
%      analyzed movie. See functions "getDiffModes" and
%      "trackDiffModeAnalysis" for detailed description of saved output.
%
%Khuloud Jaqaman, August 2017
%
% Jesus Vega-Lugo, September 2023, Updated to allow the used of clean
% tracks coming from refineTracksWithinMask. Updated documentation and use
% input parser.
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

%% Parse input

ip = inputParser;

ip.CaseSensitive = false;
ip.FunctionName = 'analyzeDiffusionModesMLMD';

addRequired(ip,'MLMD',@(x) isa(x,'MovieData') || isa(x,'MovieList'))
addParameter(ip,'channel2analyze',1,@isscalar)
addParameter(ip,'minLength',5,@isscalar)
addParameter(ip,'alpha',0.01,@isscalar)
addParameter(ip,'maxNumMode',10,@isscalar)
addParameter(ip,'binStrategy',2,@isscalar)
addParameter(ip,'subSampSize',[],@isscalar)
addParameter(ip,'doControl',1,@(x) x==1 || x==0)
addParameter(ip,'forceDecompose',1,@(x) x==1 || x==0)
addParameter(ip,'showPlot',1,@(x) x==1 || x==0)
addParameter(ip,'plotName','Figure',@(x) isstring(x) || ischar(x))
addParameter(ip,'useCleanTracks',0,@(x) x==1 || x==0)
addParameter(ip,'diffModeDividerStruct',[],@isstruct)

parse(ip,MLMD,varargin{:})

channel2analyze = ip.Results.channel2analyze;
minLength = ip.Results.minLength;
alpha = ip.Results.alpha;
maxNumMode = ip.Results.maxNumMode;
binStrategy = ip.Results.binStrategy;
subSampSize = ip.Results.subSampSize;
doControl = ip.Results.doControl;
forceDecompose = ip.Results.forceDecompose;
showPlot = ip.Results.showPlot;
plotName = ip.Results.plotName;
useCleanTracks = ip.Results.useCleanTracks;
diffModeDividerStruct = ip.Results.diffModeDividerStruct;

%% Analysis

%determine if input is a MovieList or MovieData object
if isa(MLMD,'MovieList') %if it is a movieDist
    
    listFlag = 1;
    
    %rename to ML
    ML = MLMD;
    clear MLMD
    
    %get number of movies
    numMovies = length(ML.movieDataFile_);
    
else %if it is a movieData
    
    listFlag = 0;
    
    %rename to MD
    MD = MLMD;
    clear MLMD
    
    %assign number of movies to 1
    numMovies = 1;
    
end

%go over all movies and run diffusion mode analysis
for iM = 1 : numMovies
    
    %get movieData of current movie
    if listFlag == 1
        MD = MovieData.load(ML.movieDataFile_{iM});
    end
    
    %add diffusion mode analysis process if never run before
    iProc = MD.getProcessIndex('DiffusionModeAnalysisProcess',1,0);
    if isempty(iProc)
        iProc=numel(MD.processes_)+1;
        MD.addProcess(DiffusionModeAnalysisProcess(MD));
    end
    
    %define function parameters
    
    %general
    funParams = MD.getProcess(iProc).funParams_;
    funParams.ChannelIndex = channel2analyze;
    
    %function-specific
    funParams.minLength = minLength;
    funParams.alpha = alpha;
    funParams.showPlot = showPlot;
    funParams.maxNumMode = maxNumMode;
    funParams.binStrategy = binStrategy;
    funParams.plotName = plotName;
    funParams.subSampSize = subSampSize;
    funParams.doControl = doControl;
    funParams.diffModeDividerStruct = diffModeDividerStruct;
    funParams.forceDecompose = forceDecompose;
    funParams.useCleanTracks = useCleanTracks;
    
    %general again
    parseProcessParams(MD.getProcess(iProc),funParams);
    
    %% Run analysis processes
    cellfun(@run,MD.processes_(iProc));
    
end


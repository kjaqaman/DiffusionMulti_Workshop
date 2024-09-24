function analyzeTransientDiffusionMLMD(MLMD,varargin)
%ANALYZETRANSIENTDIFFUSIONMLMD performs transient diffusion analysis on on tracks in MD/ML
%
%SYNOPSIS analyzeTransientDiffusionMLMD(MLMD,varargin)
%
%INPUT 
%       MLMD        : MovieList or MovieData object for movie(s) to be analyzed
%
%OPTIONAL (as name-value pair)
%
%       peakConfLevel: Confidence level (in percent) for choosing peaks
%                     when initially segmenting track. Default : 95.
%       probDim     : Problem dimensionality.
%                     Optional. Default: 2.
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0. 
%                     NOT RECOMMENDED IF RUNNING ON MOVIELIST.
%                     Results can be plotted only if problem is 2D.
%                     color-coding:
%                     *brown: immobile
%                     *blue: confined diffusion.
%                     *cyan: normal diffusion.
%                     *magenta: super diffusion.
%                     *black: unclassified.
%
%    useCleanTracks : 1 to use cleaned tracks, 0 to use original tracks
%                     Default: 0
%
%OUPUT 
%      Output is saved in directory TransientDiffusionAnalysis belonging to each
%      analyzed movie. See function basicTransientDiffusionAnalysisv1 for
%      detailed description of saved output.
%
%Khuloud Jaqaman, May 2019
%
% Jesus Vega-Lugo, September 2023, Updated to allow the used of clean
% tracks coming from refineTracksWithinMask. Updated to use input parser
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
ip.FunctionName = 'analyzeTransientDiffusionMLMD';

addRequired(ip,'MLMD',@(x) isa(x,'MovieData') || isa(x,'MovieList'))
addParameter(ip,'peakConfLevel',95,@isscalar)
addParameter(ip,'probDim',2,@isscalar)
addParameter(ip,'plotRes',0,@(x) x==1 || x==0)
addParameter(ip,'channel2analyze',1,@isscalar)
addParameter(ip,'useCleanTracks',0,@(x) x==1 || x==0)

parse(ip,MLMD,varargin{:});

peakConfLevel = ip.Results.peakConfLevel;
probDim = ip.Results.probDim;
plotRes = ip.Results.plotRes;
channel2analyze = ip.Results.channel2analyze;
useCleanTracks = ip.Results.useCleanTracks;


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
    
    %add transient diffusion analysis process if never run before
    iProc = MD.getProcessIndex('TransientDiffusionAnalysisProcess',1,0);
    if isempty(iProc)
        iProc=numel(MD.processes_)+1;
        MD.addProcess(TransientDiffusionAnalysisProcess(MD));
    end
    
    %define function parameters
    
    %general
    funParams = MD.getProcess(iProc).funParams_;
    funParams.ChannelIndex = channel2analyze;
    
    %function-specific
    funParams.peakAlpha = peakConfLevel;
    funParams.probDim = probDim;
    funParams.plotRes = plotRes;
    funParams.useCleanTracks = useCleanTracks;
    %general again
    parseProcessParams(MD.getProcess(iProc),funParams);
    
    %% Run analysis processes
    cellfun(@run,MD.processes_(iProc));
    
end


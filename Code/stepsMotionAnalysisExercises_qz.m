% Now using MATLAB-based testing/performance suite
% Test u-track
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

% preconditions
disp(['PWD:', pwd]);
% Dump path for debugging
s_path = strsplit(path,':');
s_match = (cellfun(@(x) regexp(x,'toolbox'), s_path, 'UniformOutput', false))';
matlab_paths = s_path(cellfun(@isempty, s_match))';
disp('    [MATLAB] current top-level paths....');
% disp(matlab_paths);

disp(['Java heap max: ' num2str(java.lang.Runtime.getRuntime.maxMemory/1e9) 'GB'])
disp('Starting u-track script');

%----Initialization of temp dir
package_name = 'utrack';
t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
mkdir(tmpdir);


cd(tmpdir);

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Test image
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zipPath = '/work/bioinformatics/s184919/Data/Khuloud/20240906RM1workshop/MotionAnalysis/VEGFR2-JF549-01.zip';

unzip(zipPath, tmpdir);

% % Analysis Output Directory - do not need this if MD created based on BioFormats
% saveFolder = [tmpdir filesep 'analysis'];
% mkdir(saveFolder);

% Imaging/Microscope Parameters
pixelSize = 81; 
timeInterval = 0.1; 
numAperture_ = 1.49;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MovieData creation

% Initialize from raw - Bioformats based
BFDataPath = [tmpdir filesep 'VEGFR2-JF549-01.tif'];
MD = MovieData(BFDataPath);

% Set some additional movie properties
MD.pixelSize_ = pixelSize;
MD.timeInterval_ = timeInterval;
MD.numAperture_ = numAperture_;

% Run sanityCheck on MovieData.
% Check image size and number of frames are consistent.
% Save the movie if successful
MD.sanityCheck;
MD.save;
MD.reset();

% Load the movie/dispaly contents
clear MD; % verify we can reload the object as intended.
MD = MovieData.load(fullfile(tmpdir,'VEGFR2-JF549-01','VEGFR2-JF549-01.mat'));


% Create InfoFlow Package and retrieve package index
Package_ = UTrackPackage(MD);
MD.addPackage(Package_);
stepNames = Package_.getProcessClassNames;
iPack =  MD.getPackageIndex('UTrackPackage');
disp('=====================================');
disp('|| Available Package Process Steps ||');
disp('=====================================');
disp(MD.getPackage(1).getProcessClassNames');

steps2Test = [1, 2, 3];
assert(length(Package_.processes_) >= length(steps2Test));
assert(length(Package_.processes_) >= max(steps2Test));
disp('Selected Package Process Steps');

for i=steps2Test
  disp(['Step ' num2str(i) ': ' stepNames{i}]);
end


%% Step 1: SubResolutionProcess
disp('===================================================================');
disp('Running (1st) Process');
disp('===================================================================');
iPack = 1;
step_ = 1;
MD.getPackage(iPack).createDefaultProcess(step_)
params = MD.getPackage(iPack).getProcess(step_).funParams_;

params.detectionParam.psfSigma = 1.3577;
params.detectionParam.alphaLocMax = 0.1;
params.detectionParam.doMMF = 1;

MD.getPackage(iPack).getProcess(step_).setPara(params);
MD.save;
params = MD.getPackage(iPack).getProcess(step_).funParams_
MD.getPackage(iPack).getProcess(step_).run();


%% Step 2: TrackingProcess
disp('===================================================================');
disp('Running (2nd) Process');
disp('===================================================================');
iPack = 1;
step_ = 2;
MD.getPackage(iPack).createDefaultProcess(step_)
params = MD.getPackage(iPack).getProcess(step_).funParams_;

params.gapCloseParam.mergeSplit = 1;
params.costMatrices(1,1).parameters.maxSearchRadius = 8;
params.costMatrices(1,1).parameters.diagnostics = [];
params.costMatrices(1,2).parameters.maxSearchRadius = 8;
params.costMatrices(1,2).parameters.useLocalDensity = 0;
params.costMatrices(1,2).parameters.brownScaling = [0.25 0.01];
params.costMatrices(1,2).parameters.resLimit = 3;
params.costMatrices(1,2).parameters.gapExcludeMS = 1;
params.costMatrices(1,2).parameters.strategyBD = -1;

MD.getPackage(iPack).getProcess(step_).setPara(params);
MD.save;
params = MD.getPackage(iPack).getProcess(step_).funParams_
MD.getPackage(iPack).getProcess(step_).run();


%% Step 3: MotionAnalysisProcess
disp('===================================================================');
disp('Running (3rd) Process');
disp('===================================================================');
iPack = 1;
step_ = 3;
MD.getPackage(iPack).createDefaultProcess(step_)
params = MD.getPackage(iPack).getProcess(step_).funParams_;

params.alphaValues = [-0.0500 0.1000];

MD.getPackage(iPack).getProcess(step_).setPara(params);
MD.save;
params = MD.getPackage(iPack).getProcess(step_).funParams_
MD.getPackage(iPack).getProcess(step_).run();


%%
% %(1) Open u-quantify GUI and visualize tracks.
% u_quantify
% 
% %(2) Run Moment Scaling Spectrum (MSS) diffusion analysis as part of u-track
% %package. This is "Step 3: Track Analysis" in the GUI.
% 
% %In the track analysis interface, select "Check to use threshold that
% %balances the error rate ...". Otherwise, use default parameters.
% 
% %(2') Visualize MSS analysis results in u-quantify GUI.
% %(2'') Open output .mat file (in uTrackPackage/MotionAnalysis subdirectory) and go over analysis summary.
% 
% %(3) Save and Exit uTrackPackageGUI.
% 
% %(4) Load the MovieData object of movie being analyzed (first navigate to
% %its subdirectory)
% MD = MovieData.load('VEGFR2-JF549-01.mat');

%(5) Run Diffusion Mode analysis
diffModeDividerStruct = struct('trajLength',5','coordStd',0.25','divider',0.3981);
analyzeDiffusionModesMLMD(MD,'forceDecompose',0, 'showPlot',0,...
    'diffModeDividerStruct',diffModeDividerStruct);

%(5') Visualize Diff Mode analysis results in u-quantify GUI.
%(5'') Open output .mat file (in DiffusionModeAnalysis subdirectory) and go over analysis summary.

%(6) Run Transient Diffusion (DC-MSS) analysis
analyzeTransientDiffusionMLMD(MD,'peakConfLevel',95,'plotRes',1);

%(6') Visualize DC-MSS results in plotted figure

%(6'') Analyze DC-MSS results

%Move all the output files into one common subdirectory called
%"transientAnalysis" (rename thee files by adding _1, _2, etc.)

%Navigate to that subdirectory

cd VEGFR2-JF549-01/TransientDiffusionAnalysis % QZ added

%Run the summary analysis
[fullLifeTotal,fullLifeTotalSE,diffProbTotal,diffRate,diffRateSE] = ...
    plotTransientStatsWrapper(pwd,1,1,1,5);

%(7) Miscellaneous visualization functions that I might use
plotTracksDiffAnalysis2D(tracks,diffAnalysisRes);
plotTracksTransDiffAnalysis2D(tracks,diffAnalysisRes);
plotTracks2D(tracks,[],'3',[],0,1);
plotCompTrack(tracks(1));




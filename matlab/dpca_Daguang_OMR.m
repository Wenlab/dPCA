%% load data and set parameters
path_working = 'E:\\WenLab\\neural_activity_analysis\\';
load(fullfile(path_working,'2025-10-21_OMR\\CalTrace\\trials.mat'));
trials = trials_control; % trials_test/trials_stimulus_control/trials_stimulus_test
trial_type = 1; % 1: trial sync with bout; 2: trial sync with stimulus

num_trials = size(trials,2);
fps_fluore = 2;
length_trial = size(trials(1,1).CalTrace,2); % (frames)
ifSimultaneousRecording = true;

%% construct firingRates and trialNum
N = size(trials(1,1).CalTrace,1);           % number of neurons
T = length_trial;                           % number of time points
S = 2;                                      % number of stimuli: left: -1->1, right: 1->2
D = 2;                                      % number of decisions: bout: 1(left: -1, right: 1), no bout: 3->2
E = num_trials;                             % maximal number of trial repetitions

if trial_type==1
    firingRates = NaN(N, S, D, T, E);
    trialNum = zeros(N,S,D);
elseif trial_type==2
    firingRates = NaN(N, S, T, E);
    trialNum = zeros(N,S);
end

for i=1:num_trials
    if trials(1,i).Stimulus==-1
        ss = 1;
    elseif trials(1,i).Stimulus==1
        ss = 2;
    end
    if trial_type==1
        if abs(trials(1,i).Decision)==1
            dd = 1;
        elseif trials(1,i).Decision==3
            dd = 2;
        end
        trialNum(:,ss,dd) = trialNum(:,ss,dd)+1;
        ee = trialNum(1,ss,dd);
        firingRates(:,ss,dd,:,ee) = reshape(trials(1,i).CalTrace,N,1,1,T,1);
    elseif trial_type==2
        trialNum(:,ss) = trialNum(:,ss)+1;
        ee = trialNum(1,ss);
        firingRates(:,ss,:,ee) = reshape(trials(1,i).CalTrace,N,1,T,1);
    end
end

E = max(trialNum,[],'all'); % correct E
if trial_type==1
    firingRates = firingRates(:,:,:,:,1:E);
elseif trial_type==2
    firingRates = firingRates(:,:,:,1:E);
end

% delete the stimulus or decision type with less than 2 trial

if trial_type==1
    firingRatesAverage = mean(firingRates, 5,'omitnan');
elseif trial_type==2
    firingRatesAverage = mean(firingRates, 4,'omitnan');
end

%% Define parameter grouping

% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; here we ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:
if trial_type==1
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};
elseif trial_type==2
    combinedParams = {{1, [1 2]}, {2}};
    margNames = {'Stimulus', 'Condition-independent'};
    margColours = [23 100 171; 150 150 150]/256;
end

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
% timeEvents = time(round(length(time)/2));

% check consistency between trialNum and firingRates
if trial_type==1
    for n = 1:size(firingRates,1)
        for s = 1:size(firingRates,2)
            for d = 1:size(firingRates,3)
                assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
            end
        end
    end
elseif trial_type==2
    for n = 1:size(firingRates,1)
        for s = 1:size(firingRates,2)
            assert(isempty(find(isnan(firingRates(n,s,:,1:trialNum(n,s))), 1)), 'Something is wrong!')
        end
    end
end

%% dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'noiseCovType', 'averaged', ...
    'lambdas', 1e-03 * 3.23 * 1.5.^[0:22], ...
    'numRep', 10, ...  % increase this number to ~10 for better accuracy
    'filename', fullfile(path_working,'tmp_optimalLambdas.mat'));

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording,'type','averaged');

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);


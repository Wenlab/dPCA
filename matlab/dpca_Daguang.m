%% load data and set parameters
path_working = 'D:\\Nutstore\\�ҵļ����\\��ʱ\\2023_11_28-16_56_8\\';
load(fullfile(path_working,'neural_activity\\CalTrace\\prebout_activity.mat'));

num_bouts = length(prebout_activity);

fps_fluore = 2;
length_pre = 1;
length_post = 1; 
length_trial = length_pre + length_post + 1; % (frames)
ifSimultaneousRecording = true;

%% extract neural activity of each bout

N = size(prebout_activity(1).activity,1);   % number of neurons
T = length_trial;                           % number of time points
S = 1;                                      % number of stimuli: spontaneous
D = 13;                                     % number of decisions: 13 types of bout
E = num_bouts;                              % maximal number of trial repetitions

% firingRates = NaN(N, S, D, T, E);
% trialNum = zeros(N,S,D);
firingRates = NaN(N, D, T, E);
trialNum = zeros(N,D);

for i=1:num_bouts
    % ss = 1;
    dd = abs(prebout_activity(i).bout_type);
    % trialNum(:,ss,dd) = trialNum(:,ss,dd)+1;
    % ee = trialNum(1,ss,dd);
    % firingRates(:,ss,dd,:,ee) = reshape(prebout_activity(i).activity,N,1,1,T,1);
    trialNum(:,dd) = trialNum(:,dd)+1;
    ee = trialNum(1,dd);
    activity_trial = Cal_spike_smooth(:,[prebout_activity(i).frame_index-length_pre:prebout_activity(i).frame_index+length_post]);
    if isempty(find(isnan(activity_trial), 1))
        firingRates(:,dd,:,ee) = reshape(activity_trial,N,1,T,1);
    else
        trialNum(:,dd) = trialNum(:,dd)-1;
        continue;
    end
end

E = max(trialNum,[],'all'); % correct E
% firingRates = firingRates(:,:,:,:,1:E);
firingRates = firingRates(:,:,:,1:E);

% delete the bout type with less than 2 trial
firingRates(:,[3,4,10],:,:) = [];
trialNum(:,[3,4,10]) = [];

% firingRatesAverage = mean(firingRates, 5,'omitnan');
firingRatesAverage = mean(firingRates, 4,'omitnan');

%% Define parameter grouping

% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
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

% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
combinedParams = {{1, [1 2]}, {2}};
margNames = {'Decision', 'Condition-independent'};
margColours = [187 20 25; 150 150 150]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
% timeEvents = time(round(length(time)/2));

% check consistency between trialNum and firingRates
% for n = 1:size(firingRates,1)
%     for s = 1:size(firingRates,2)
%         for d = 1:size(firingRates,3)
%             assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
%         end
%     end
% end
for n = 1:size(firingRates,1)
    for d = 1:size(firingRates,2)
        assert(isempty(find(isnan(firingRates(n,d,:,1:trialNum(n,d))), 1)), 'Something is wrong!')
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


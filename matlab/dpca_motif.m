%% load data and set parameters
path_working = 'D:\\Nutstore\\我的坚果云\\临时\\2023_11_28-16_56_8\\';
load(fullfile(path_working,'neural_activity\\CalTrace\\prebout_activity.mat'));
load(fullfile(path_working,'behavior\\bout_seq_segmentation_and_motifs_by_interval.mat'));

motif_picked = [4,5,6,8,10,20,22,23,24];
fps_fluore = 2;
length_trial = 2; % (frames) (which is equal to the length of these motifs)
ifSimultaneousRecording = true;

%% extract neural activity of each bout

N = size(prebout_activity(1).activity,1);                   % number of neurons
T = length_trial;                                           % number of time points
S = 1;                                                      % number of stimuli: spontaneous
D = length(motif_picked);                                   % number of decisions: num of different motifs
E = max(cat(1,motifs_sort(motif_picked).num_occurrences));  % maximal number of trial repetitions

% firingRates = NaN(N, S, D, T, E);
% trialNum = zeros(N,S,D);
firingRates = NaN(N, D, T, E);
trialNum = zeros(N,D);

for dd=1:D
    motif_present = motif_picked(dd);
    seg_idx_motif_present = motifs_sort(motif_present).seg_idx;
    for ee=1:length(seg_idx_motif_present)
        bout_time_fluo_seg_present = seg_by_interval(seg_idx_motif_present(ee)).bout_time_fluo;
        activity_trial = Cal_spike_smooth(:,bout_time_fluo_seg_present);
        if isempty(find(isnan(activity_trial), 1))
            firingRates(:,dd,:,ee) = reshape(activity_trial,N,1,T,1);
            trialNum(:,dd) = trialNum(:,dd)+1;
        else
            continue;
        end
    end
end

E = max(trialNum,[],'all'); % correct E
% firingRates = firingRates(:,:,:,:,1:E);
firingRates = firingRates(:,:,:,1:E);

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
margNames = {'Motif', 'Condition-independent'};
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

% optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
%     'numComps', T*S*D-1, ...
%     'combinedParams', combinedParams, ...
%     'simultaneous', ifSimultaneousRecording, ...
%     'noiseCovType', 'averaged', ...
%     'lambdas', 1e-03 * 3.23 * 1.5.^[0:22], ...
%     'numRep', 10, ...  % increase this number to ~10 for better accuracy
%     'filename', fullfile(path_working,'tmp_optimalLambdas.mat'));

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording,'type','averaged');

lambda = 0.01;
[W,V,whichMarg] = dpca(firingRatesAverage, T*S*D-1, ...
    'combinedParams', combinedParams, ...
    'lambda', lambda, ...
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


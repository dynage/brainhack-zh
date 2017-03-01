%% Tutorial Code for Brainhack Zurich 2017, Microstates EEGLAB TOOLBOX
clear;
addpath('~/Dropbox/EEG_analysis/GeneralMatlab/eeglab14_0_0b');
eeglab;
close;
%% load the EEG file for the segmentation
% this is resting EEG filtered at 2 - 20 Hz
load('~/dropbox/Brainhack/PreprocessedData/gp_m10916141001.mat');

% force average reference
mEEG = mean(EEG.data,1);
EEG.data = EEG.data - repmat(mEEG,size(EEG.data,1),1);

%% Parse GFP peaks
% we do the segmentation on the GFP peaks, because they are high in
% signal-to-noise

% calculate the GFP
GFP = double(std(EEG.data,[],1));
[~, peakidx] = findpeaks(GFP,'MinPeakDistance',15); % 15 = minimum distance between two peaks
% check if the peak searching algorithm does something appropriate
findpeaks(GFP(1:500),'MinPeakDistance',15,'Annotate','extents');
% select a number of random peaks (if there are many peaks or you are
% impatient)
selection = randperm(length(peakidx));
EEG.data = EEG.data(:,peakidx(selection(1:1000))); % here I select 1000 peaks

%% Microstate Segmentation
% we use modified K-means with 4 microstates, by default modified k-means
% ignored polarity. See help pop_micro_segment for more information about
% defaults and options
OUTEEG = pop_micro_segment(EEG,'algorithm','Modified K-means','Nmicrostates',4);

% store the prototype topographies into a new variable, which will serve as
% template for the backfitting
MStemp = OUTEEG.microstate.scalp_maps;

%% create figure with template maps
% we want to see how the template maps look like
f = figure
k = 4
for i = 1:k
subplot(1,k,i)   ;
topoplot(MStemp(:,i),EEG.chanlocs,'electrodes','off','style','map','headrad',0);
end

%% We reload the same EEG file for the backfitting
% here we use the whole file (not only the GFP peaks)
load('~/dropbox/Brainhack/PreprocessedData/gp_m10916141001.mat');

%% epoch data
%EEG = pop_epoch(EEG, {20,30}, [0 20]);

%% fit back to the EEG data
% we use MStemp as prototype maps choose that a microstate should be not
% less than 15 TF (=30 ms), we ignore polarity and do not select the option
% to sequentialize
% The Output will create the structure EEG.microstate.fit
[EEG] = MicroFit(EEG,MStemp,'minTF',15,'polarity',0,'sequentialize',1);

%% Statistics
% here we calculate the statistics of the backfitting
% The output will create the structure EEG.microstate.stats
[EEG] =  MicroStats(EEG);

%% Plot the segmented GFP
MicroPlot(EEG,'epoch',1:10000);
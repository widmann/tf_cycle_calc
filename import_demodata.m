function TMP = import_demodata()
% import_demodata - Download and convert Cohen, 2019 demo data to EEGLAB format
% Copyright (c) 2025 Andreas Widmann, University of Leipzig
% Author: Andreas Widmann, widmann@uni-leipzig.de
%
% Refs:
%   [1] Cohen, M.X. (2019). A better way to define and describe Morlet
%       wavelets for time-frequency analysis, NeuroImage, 199, 81-86.

% Download and load Cohen, 2019 demo data
if exist('demodata_cohen2019.set', 'file')

    TMP = pop_loadset('filename', 'demodata_cohen2019.set');

else

    if ~exist('MorletWaveletDefinition_data.mat', 'file')
        url = 'https://ars.els-cdn.com/content/image/1-s2.0-S1053811919304409-mmc2.mat';
        filename = 'MorletWaveletDefinition_data.mat';
        websave(filename, url);
    end
    load('MorletWaveletDefinition_data.mat', 'EEG')

    % Convert to EEGLAB format
    TMP = eeg_emptyset;

    TMP.data(1, :, :) = EEG.data;
    TMP.nbchan = size(TMP.data, 1);
    TMP.pnts = size(TMP.data, 2);
    TMP.trials = size(TMP.data, 3);

    TMP.srate = EEG.srate;
    TMP.xmin = EEG.times(1) / 1000;
    TMP.xmax = (TMP.pnts - 1) / TMP.srate + TMP.xmin;

    TMP = eeg_checkset(TMP);
    pop_saveset(TMP, 'filename', 'demodata_cohen2019.set');

end

end
function demo_cohen2019_eeglab()
% demo_cohen2019_eeglab - Demonstrate tf_cycle_calc with Cohen, 2019 demo
% data and EEGLAB pop_newtimef
% Copyright (c) 2025 Andreas Widmann, University of Leipzig
% Author: Andreas Widmann, widmann@uni-leipzig.de
%
% Refs:
%   [1] Cohen, M.X. (2019). A better way to define and describe Morlet
%       wavelets for time-frequency analysis, NeuroImage, 199, 81-86.

EEG = import_demodata;

freqs = 4:40;
basewin = [-0.5 -0.2];

hFig = figure;


%%%% STFT with 0.1 s Gaussian window (Gabor transform, STFT)

fwhm_t = 0.1; % Constant duration
cycles = tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

figure(hFig), hPlot(1) = subplot(2, 3, 1);
pop_newtimef(EEG, 1, 1, [-650 1350], cycles, 'baseline', basewin * 1000, 'freqs', freqs, 'erspmax', [-2.5 2.5], 'plotitc' , 'off', 'plotphase', 'off', 'ntimesout', 400);
title(hPlot(1), 'A: Constant FWHM width 0.1 s (STFT)', 'FontSize', 10)


%%%% STFT with 0.5 s Gaussian window (Gabor transform, STFT)

fwhm_t = 0.5; % Constant duration
cycles = tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

figure(hFig), hPlot(2) = subplot(2, 3, 2);
pop_newtimef(EEG, 1, 1, [-1244 1944], cycles, 'baseline', basewin * 1000, 'freqs', freqs, 'erspmax', [-2.5 2.5], 'plotitc' , 'off', 'plotphase', 'off', 'ntimesout', 400);
title(hPlot(2), 'B: Constant FWHM width 0.5 s (Gabor)', 'FontSize', 10)


%%%% Fixed 7 cycles

cycles = 7;
cycles = tf_cycle_calc('freqs', freqs, 'width', cycles, 'width_unit', 'cycles', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

figure(hFig), hPlot(3) = subplot(2, 3, 3);
pop_newtimef(EEG, 1, 1, [-1475 2175], cycles, 'baseline', basewin * 1000, 'freqs', freqs, 'erspmax', [-2.5 2.5], 'plotitc' , 'off', 'plotphase', 'off', 'ntimesout', 400);
title(hPlot(3), 'C: 7 cycles (Morlet wavelet)', 'FontSize', 10)


%%%% Variable cycles

cycles = linspace(3, 8, length(freqs));
cycles = tf_cycle_calc('freqs', freqs, 'width', cycles, 'width_unit', 'cycles', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

figure(hFig), hPlot(4) = subplot(2, 3, 4);
pop_newtimef(EEG, 1, 1, [-918 1618], cycles, 'baseline', basewin * 1000, 'freqs', freqs, 'erspmax', [-2.5 2.5], 'plotitc' , 'off', 'plotphase', 'off', 'ntimesout', 400);
title(hPlot(4), 'D: 3 to 8 cycles (linear)', 'FontSize', 10)


%%%% Variable temporal width FWHM from 0.5 to 0.2 s

fwhm_t = linspace(0.5, 0.2, length(freqs));
cycles = tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

figure(hFig), hPlot(5) = subplot(2, 3, 5);
pop_newtimef(EEG, 1, 1, [-1244 1944], cycles, 'baseline', basewin * 1000, 'freqs', freqs, 'erspmax', [-2.5 2.5], 'plotitc' , 'off', 'plotphase', 'off', 'ntimesout', 400);
title(hPlot(5), 'E: 0.5 to 0.2 s FWHM width (linear)', 'FontSize', 10)


%%%% Variable spectral bandwidth FWHM from 2 to 6 Hz

fwhm_f = linspace(2, 6, length(freqs));
cycles = tf_cycle_calc('freqs', freqs, 'width', fwhm_f, 'width_unit', 'fwhm_f', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

figure(hFig), hPlot(6) = subplot(2, 3, 6);
pop_newtimef(EEG, 1, 1, [-1155 1855], cycles, 'baseline', basewin * 1000, 'freqs', freqs, 'erspmax', [-2.5 2.5], 'plotitc' , 'off', 'plotphase', 'off', 'ntimesout', 400);
title(hPlot(6), 'F: 2 to 6 Hz FWHM bandwidth (linear)', 'FontSize', 10)


% Scaling, labels, export, etc.
set(findobj('Tag', 'AxCycleCalcCycles'), 'YLim', [0 40])
h = findobj('Tag', 'AxCycleCalcWidth');
legend(h, get([hPlot.Title], 'String'))

exportgraphics(hFig, 'demo_cohen2019_eeglab-1.png', 'Resolution', '75')
exportgraphics(findobj('tag', 'FigCycleCalc'), 'demo_cohen2019_eeglab-2.png', 'Resolution', '75')

end

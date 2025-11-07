function demo_cohen2019()
% demo_cohen2019 - Demonstrate tf_cycle_calc with Cohen, 2019 demo data
% Copyright (c) 2025 Andreas Widmann, University of Leipzig
% Author: Andreas Widmann, widmann@uni-leipzig.de
%
% Refs:
%   [1] Cohen, M.X. (2019). A better way to define and describe Morlet
%       wavelets for time-frequency analysis, NeuroImage, 199, 81-86.

EEG = import_demodata;

freqs = 2:40;
basewin = [-0.5 -0.2];
basewin = round((basewin - EEG.xmin) * EEG.srate + 1);
fCtr = [10 21 32];
fCtrIdx = find(ismember(freqs, fCtr));

hFig = figure;
colmap = my_brewerburd;
% colmap = turbo;
colormap(colmap)
yscale = 'linear';
% yscale = 'log';


%%%% STFT with 0.1 s Gaussian window (Gabor transform, STFT)

fwhm_t = 0.1; % Constant duration
tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

for iTrial = 1:EEG.trials
    [cycles, widths_table, tf_data(:, :, iTrial)] = tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 0, 'data', EEG.data(1, :, iTrial), 'fs', EEG.srate, 'norm', 0);
end

% Average over trials and baseline correct
tf_data = mean(abs(tf_data).^2, 3);
tf_data = 10 * log10(tf_data ./ mean(tf_data(basewin(1):basewin(2), :), 1));

figure(hFig), hPlot(1) = subplot(2, 3, 1);
contourf(EEG.times / 1000, freqs, tf_data', 64, 'LineColor', 'none')
hold('on')
for idx = 1:length(fCtrIdx)
    patch('XData', [0 widths_table(fCtrIdx(idx), 4) widths_table(fCtrIdx(idx), 4) 0], 'YData', [fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2], 'FaceColor', 'none')
end
title('A: Constant FWHM width 0.1 s (STFT)')


%%%% STFT with 0.5 s Gaussian window (Gabor transform, STFT)

fwhm_t = 0.5; % Constant duration
tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

for iTrial = 1:EEG.trials
    [cycles, widths_table, tf_data(:, :, iTrial)] = tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 0, 'data', EEG.data(1, :, iTrial), 'fs', EEG.srate, 'norm', 0);
end

% Average over trials and baseline correct
tf_data = mean(abs(tf_data).^2, 3);
tf_data = 10 * log10(tf_data ./ mean(tf_data(basewin(1):basewin(2), :), 1));

figure(hFig), hPlot(2) = subplot(2, 3, 2);
contourf(EEG.times / 1000, freqs, tf_data', 64, 'LineColor', 'none')
hold('on')
for idx = 1:length(fCtrIdx)
    patch('XData', [0 widths_table(fCtrIdx(idx), 4) widths_table(fCtrIdx(idx), 4) 0], 'YData', [fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2], 'FaceColor', 'none')
end
title('B: Constant FWHM width 0.5 s (STFT)')


%%%% Fixed 7 cycles

cycles = 7;
tf_cycle_calc('freqs', freqs, 'width', cycles, 'width_unit', 'cycles', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

for iTrial = 1:EEG.trials
    [cycles, widths_table, tf_data(:, :, iTrial)] = tf_cycle_calc('freqs', freqs, 'width', cycles, 'width_unit', 'cycles', 'plot', 0, 'data', EEG.data(1, :, iTrial), 'fs', EEG.srate, 'norm', 0);
end

% Average over trials and baseline correct
tf_data = mean(abs(tf_data).^2, 3);
tf_data = 10 * log10(tf_data ./ mean(tf_data(basewin(1):basewin(2), :), 1));

figure(hFig), hPlot(3) = subplot(2, 3, 3);
contourf(EEG.times / 1000, freqs, tf_data', 64, 'LineColor', 'none')
hold('on')
for idx = 1:length(fCtrIdx)
    patch('XData', [0 widths_table(fCtrIdx(idx), 4) widths_table(fCtrIdx(idx), 4) 0], 'YData', [fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2], 'FaceColor', 'none')
end
title('C: 7 cycles (Morlet wavelet)')


%%%% Variable cycles

cycles = linspace(3, 8, length(freqs));
tf_cycle_calc('freqs', freqs, 'width', cycles, 'width_unit', 'cycles', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

for iTrial = 1:EEG.trials
    [cycles, widths_table, tf_data(:, :, iTrial)] = tf_cycle_calc('freqs', freqs, 'width', cycles, 'width_unit', 'cycles', 'plot', 0, 'data', EEG.data(1, :, iTrial), 'fs', EEG.srate, 'norm', 0);
end

% Average over trials and baseline correct
tf_data = mean(abs(tf_data).^2, 3);
tf_data = 10 * log10(tf_data ./ mean(tf_data(basewin(1):basewin(2), :), 1));

figure(hFig), hPlot(4) = subplot(2, 3, 4);
contourf(EEG.times / 1000, freqs, tf_data', 64, 'LineColor', 'none')
hold('on')
for idx = 1:length(fCtrIdx)
    patch('XData', [0 widths_table(fCtrIdx(idx), 4) widths_table(fCtrIdx(idx), 4) 0], 'YData', [fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2], 'FaceColor', 'none')
end
title('D: 3 to 8 cycles (linear)')


%%%% Variable temporal width FWHM from 0.5 to 0.2 s

fwhm_t = linspace(0.5, 0.2, length(freqs));
tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

for iTrial = 1:EEG.trials
    [cycles, widths_table, tf_data(:, :, iTrial)] = tf_cycle_calc('freqs', freqs, 'width', fwhm_t, 'width_unit', 'fwhm_t', 'plot', 0, 'data', EEG.data(1, :, iTrial), 'fs', EEG.srate, 'norm', 0);
end

% Average over trials and baseline correct
tf_data = mean(abs(tf_data).^2, 3);
tf_data = 10 * log10(tf_data ./ mean(tf_data(basewin(1):basewin(2), :), 1));

figure(hFig), hPlot(5) = subplot(2, 3, 5);
contourf(EEG.times / 1000, freqs, tf_data', 64, 'LineColor', 'none')
hold('on')
for idx = 1:length(fCtrIdx)
    patch('XData', [0 widths_table(fCtrIdx(idx), 4) widths_table(fCtrIdx(idx), 4) 0], 'YData', [fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2], 'FaceColor', 'none')
end
title('E: 0.5 to 0.2 s FWHM width (linear)')


%%%% Variable spectral bandwidth FWHM from 2 to 6 Hz

fwhm_f = linspace(2, 6, length(freqs));
tf_cycle_calc('freqs', freqs, 'width', fwhm_f, 'width_unit', 'fwhm_f', 'plot', 1); % Compare spectral bandwidth, temporal width, cycles

for iTrial = 1:EEG.trials
    [cycles, widths_table, tf_data(:, :, iTrial), cmorlf] = tf_cycle_calc('freqs', freqs, 'width', fwhm_f, 'width_unit', 'fwhm_f', 'plot', 0, 'data', EEG.data(1, :, iTrial), 'fs', EEG.srate, 'norm', 0);
end

% Average over trials and baseline correct
tf_data = mean(abs(tf_data).^2, 3);
tf_data = 10 * log10(tf_data ./ mean(tf_data(basewin(1):basewin(2), :), 1));

figure(hFig), hPlot(6) = subplot(2, 3, 6);
contourf(EEG.times / 1000, freqs, tf_data', 64, 'LineColor', 'none')
hold('on')
for idx = 1:length(fCtrIdx)
    patch('XData', [0 widths_table(fCtrIdx(idx), 4) widths_table(fCtrIdx(idx), 4) 0], 'YData', [fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) - widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2 fCtr(idx) + widths_table(fCtrIdx(idx), 3) / 2], 'FaceColor', 'none')
end
title('F: 2 to 6 Hz FWHM bandwidth (linear)')


% Scaling, labels, export, etc.
set(hPlot, 'XLim', [-0.5 1.3], 'CLim', [-2.5 2.5], 'YScale', yscale)
set([hPlot.XLabel], 'String', 'Time [s]')
set([hPlot.YLabel], 'String', 'Frequency [Hz]')
for iAx = 1:6, hCol = colorbar(hPlot(iAx)); set([hCol.Label], 'String', 'Power change [dB]'); end

set(findobj('Tag', 'AxCycleCalcCycles'), 'YLim', [0 40])
h = findobj('Tag', 'AxCycleCalcWidth');
legend(h, get([hPlot.Title], 'String'))

exportgraphics(hFig, 'demo_cohen2019-1.png', 'Resolution', '75')
exportgraphics(findobj('tag', 'FigCycleCalc'), 'demo_cohen2019-2.png', 'Resolution', '75')

% figure
% plot(cmorlf)

end

function colmap = my_brewerburd
% Refs:
%   https://colorbrewer2.org/#type=diverging&scheme=RdBu
%   https://github.com/DrosteEffect/BrewerMap
%
% Note:
%   Included here to reduce dependencies

colmap = [
    0.0196    0.1882    0.3804; ...
    0.0355    0.2242    0.4418; ...
    0.0520    0.2594    0.4996; ...
    0.0691    0.2937    0.5523; ...
    0.0868    0.3270    0.5987; ...
    0.1050    0.3594    0.6375; ...
    0.1237    0.3908    0.6675; ...
    0.1425    0.4209    0.6887; ...
    0.1607    0.4494    0.7053; ...
    0.1792    0.4769    0.7191; ...
    0.1989    0.5039    0.7312; ...
    0.2208    0.5311    0.7431; ...
    0.2458    0.5591    0.7559; ...
    0.2757    0.5885    0.7710; ...
    0.3163    0.6203    0.7879; ...
    0.3659    0.6535    0.8057; ...
    0.4205    0.6869    0.8237; ...
    0.4764    0.7191    0.8412; ...
    0.5298    0.7489    0.8575; ...
    0.5769    0.7750    0.8719; ...
    0.6205    0.7986    0.8854; ...
    0.6636    0.8211    0.8985; ...
    0.7052    0.8421    0.9109; ...
    0.7446    0.8616    0.9222; ...
    0.7810    0.8793    0.9320; ...
    0.8136    0.8951    0.9399; ...
    0.8428    0.9093    0.9461; ...
    0.8704    0.9227    0.9520; ...
    0.8963    0.9351    0.9574; ...
    0.9201    0.9465    0.9620; ...
    0.9415    0.9565    0.9656; ...
    0.9603    0.9650    0.9680; ...
    0.9726    0.9635    0.9567; ...
    0.9793    0.9510    0.9310; ...
    0.9844    0.9357    0.9031; ...
    0.9882    0.9180    0.8732; ...
    0.9906    0.8982    0.8415; ...
    0.9919    0.8769    0.8081; ...
    0.9921    0.8541    0.7732; ...
    0.9904    0.8276    0.7339; ...
    0.9866    0.7969    0.6903; ...
    0.9810    0.7628    0.6444; ...
    0.9741    0.7264    0.5983; ...
    0.9662    0.6888    0.5541; ...
    0.9577    0.6508    0.5136; ...
    0.9469    0.6120    0.4765; ...
    0.9314    0.5707    0.4402; ...
    0.9124    0.5276    0.4053; ...
    0.8913    0.4834    0.3721; ...
    0.8692    0.4386    0.3410; ...
    0.8475    0.3941    0.3125; ...
    0.8272    0.3489    0.2865; ...
    0.8076    0.2988    0.2609; ...
    0.7877    0.2467    0.2364; ...
    0.7668    0.1960    0.2138; ...
    0.7440    0.1503    0.1939; ...
    0.7182    0.1132    0.1775; ...
    0.6887    0.0874    0.1653; ...
    0.6534    0.0654    0.1545; ...
    0.6123    0.0449    0.1445; ...
    0.5662    0.0269    0.1357; ...
    0.5157    0.0127    0.1286; ...
    0.4614    0.0034    0.1237; ...
    0.4039         0    0.1216 ...
];

end
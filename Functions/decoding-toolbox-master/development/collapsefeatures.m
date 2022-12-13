clear all;

addpath('~/-pimmos/decoding/');

% Gives me: X (labels), Y (data), time
tmp = load('../examples/testdata_2class.mat');

Y = tmp.Y;              % Y: features x time x trials
X = tmp.X';
time = tmp.time;
label = tmp.label;

numF = size(Y, 1);
numT = size(Y, 2);
numN = size(Y, 3);

% Select time points, sensors, and collapse features
if (0)          % Manually
    selChan = find(ismember(label, ft_channelselection({'M*O*', 'M*P*'}, label)));
    selTime = find(inrange(time, [.06, .3]));

    fprintf('chan: %g, time: %g, memory: %.2f gb\n', length(selChan), length(selTime), (length(selChan)*length(selTime))^2*8/1024^3);
else            % Based on separate within-time and within-sensor decoding
    numF = 7500;
    
    tmp_spat = load('pCorrect_spat_2class.mat');
    tmp_time = load('pCorrect_time_2class.mat');
    
    [~, i_spat] = sort(tmp_spat.pCorrect);
    [~, i_time] = sort(tmp_time.pCorrect);
    
    selChan = i_spat((end-round(sqrt(numF))+1):end);
    selTime = i_time((end-round(sqrt(numF))+1):end);
    
    fprintf('chan: %g, time: %g, memory: %.2f gb\n', length(selChan), length(selTime), (length(selChan)*length(selTime))^2*8/1024^3);
end    
    
Ycol = reshape(Y(selChan, selTime, :), [length(selChan)*length(selTime), numN]);

% Calculate gamma
gamma = shrinkageGamma(Ycol, 1);
%gamma = 0.02;

% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

% Decode using cross-validation
cfg = [];
cfg.feedback = 'yes';
cfg.folds = folds;
cfg.trainfun = 'train_LDA';
cfg.traincfg.gamma = gamma;
cfg.decodefun = 'decode_LDA';

Xhat = decodeCrossValidation(cfg, X, Ycol);

% Plot
m0 = mean(Xhat(X==0));
m1 = mean(Xhat(X==1));

% Classify
class = (Xhat > 0);

pCorrect = mean(X == class);

% Test and plot
[~, p] = ttest2(Xhat(X==0), Xhat(X==1));

edges = linspace(min(Xhat), max(Xhat), 15);
n0 = histc(Xhat(X==0), edges);
n1 = histc(Xhat(X==1), edges);

figure; hold on; grid on;
plot(edges, n0, 'black');
plot(edges, n1, 'red');
vabline(m0, ':black');
vabline(m1, ':red');
title(sprintf('mean diff: %.2f, log10(p): %.2f\npCorrect: %.2f', m1 - m0, log10(p), pCorrect));


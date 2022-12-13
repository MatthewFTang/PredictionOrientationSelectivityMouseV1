clear all;

addpath('../');

% Gives me: X (labels), Y (data), time
tmp = load('../examples/testdata_2class.mat');

Y = tmp.Y;              % Y: features x time x trials
X = tmp.X';
time = tmp.time;
label = tmp.label;

numF = size(Y, 1);
numT = size(Y, 2);
numN = size(Y, 3);

%% 1. No time dimension
sel_t = find(time >= .12, 1);

%% 1.1. Split into train- and test-set
numN_train = 260;
numN_test = numN - numN_train;

i0 = find(X==0); i0 = i0(randperm(length(i0), numN_train/2));
i1 = find(X==1); i1 = i1(randperm(length(i1), numN_train/2));

i_train = [i0, i1];
i_test = setdiff(1:numN, i_train)';

X_train = X(i_train);
Y_train = squeeze(Y(:, sel_t, i_train));

X_test = X(i_test);
Y_test = squeeze(Y(:, sel_t, i_test));

% Create design matrix
design = designMatrix_dummy([], X_train);

% Train
cfg = [];
cfg.gamma = 0.05;
cfg.demean = 'no';
decoderBF = train_beamformer(cfg, design, Y_train);
decoderLDA = train_LDA(cfg, X_train, Y_train);

% Decode continuous output
cfg = [];
XhatBF = decode_beamformer(cfg, decoderBF, Y_test);
XhatLDA = decode_LDA(cfg, decoderLDA, Y_test);

% Classify
classBF = (XhatBF(1, :) < XhatBF(2, :));
classLDA = XhatLDA > 0;

pCorrectBF = mean(X_test == classBF)
pCorrectLDA = mean(X_test == classLDA)

%% 1.2. Cross-validation
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

% Decode using cross-validation
cfg = [];
cfg.feedback = 'yes';
cfg.folds = folds;
cfg.trainfun = 'train_beamformer';
cfg.traincfg.gamma = 0.05;
cfg.traincfg.demean = 'no';
cfg.decodefun = 'decode_beamformer';

design = designMatrix_dummy([], X);

Xhat = decodeCrossValidation(cfg, design, squeeze(Y(:, sel_t, :)));

m0 = mean(Xhat(:, X==0), 2);
m1 = mean(Xhat(:, X==1), 2);

class = (Xhat(1, :) < Xhat(2, :));

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
title(sprintf('mean diff: %.2f, m1log10(p): %.2f\npCorrect: %.2f', m1 - m0, log10(p), pCorrect));

%% 2. Including temporal dimension
%% 2.1. Matching testing and decoding time
%% 2.1.1. Split into train- and test-set
numN_train = 260;
numN_test = numN - numN_train;

i0 = find(X==0); i0 = i0(randperm(length(i0), numN_train/2));
i1 = find(X==1); i1 = i1(randperm(length(i1), numN_train/2));

i_train = [i0, i1];
i_test = setdiff(1:numN, i_train)';

X_train = X(i_train);
Y_train = squeeze(Y(:, :, i_train));

X_test = X(i_test);
Y_test = squeeze(Y(:, :, i_test));
    
% Train
cfg = [];
cfg.feedback = 'yes';
cfg.trainfun = 'train_LDA';
cfg.traincfg.gamma = 0.05;
decoder = train_array(cfg, X_train, Y_train);

% Decode continuous output
cfg = [];
cfg.feedback = 'yes';
cfg.decodefun = 'decode_LDA';
Xhat = decode_array(cfg, decoder, Y_test);

% Classify
class = (Xhat > 0);

% Plot
m0 = squeeze(mean(Xhat(:, :, X_test==0), 3));
m1 = squeeze(mean(Xhat(:, :, X_test==1), 3));

[~, p] = ttest2(Xhat(:, :, X_test==0), Xhat(:, :, X_test==1), 'dim', 3);
p = squeeze(p);

pCorrect = squeeze(mean(class == repmat(permute(X_test, [1, 3, 2]), [1, numT, 1]), 3));

figure;
subplot(3, 1, 1); hold on; grid on;
plot(time, m0, 'black');
plot(time, m1, 'red');
subplot(3, 1, 2); 
plot(time, log10(p)); grid on;
subplot(3, 1, 3); 
plot(time, pCorrect); grid on;

%% 2.1.2. Decode using cross-validation
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

cfg = [];
cfg.folds = folds;
cfg.feedback = 'yes';
cfg.trainfun = 'train_array';
cfg.traincfg.feedback = 'yes';
cfg.traincfg.trainfun = 'train_LDA';
cfg.traincfg.traincfg.gamma = 0.05;
cfg.decodefun = 'decode_array';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_LDA';

Xhat = decodeCrossValidation(cfg, X, Y);
class = (Xhat > 0);

% Plot
m0 = squeeze(mean(Xhat(:, :, X==0), 3));
m1 = squeeze(mean(Xhat(:, :, X==1), 3));

[~, p] = ttest2(Xhat(:, :, X==0), Xhat(:, :, X==1), 'dim', 3);
p = squeeze(p);

pCorrect = squeeze(mean(class == repmat(permute(X, [1, 3, 2]), [1, numT, 1]), 3));

figure;
subplot(3, 1, 1); hold on; grid on;
plot(time, m0, 'black');
plot(time, m1, 'red');
subplot(3, 1, 2); 
plot(time, log10(p)); grid on;
subplot(3, 1, 3); 
plot(time, pCorrect); grid on;

if (0)
    save('development/pCorrect_time_2class.mat', 'pCorrect');
end

%% 2.2. Temporal generalization
%% 2.2.1. Split into train- and test-set
numN_train = 260;
numN_test = numN - numN_train;

i0 = find(X==0); i0 = i0(randperm(length(i0), numN_train/2));
i1 = find(X==1); i1 = i1(randperm(length(i1), numN_train/2));

i_train = [i0, i1];
i_test = setdiff(1:numN, i_train)';

X_train = X(i_train);
Y_train = squeeze(Y(:, :, i_train));

X_test = X(i_test);
Y_test = squeeze(Y(:, :, i_test));
    
% Train
cfg = [];
cfg.feedback = 'yes';
cfg.trainfun = 'train_LDA';
cfg.traincfg.gamma = 0.05;
decoder = train_array(cfg, X_train, Y_train);

% Decode continuous output
cfg = [];
cfg.feedback = 'yes';
cfg.decodefun = 'decode_LDA';
Xhat = decode_arrayGeneralization(cfg, decoder, Y_test);

% Classify
class = (Xhat > 0);

% Plot
m0 = squeeze(mean(Xhat(:, :, :, X_test==0), 4));
m1 = squeeze(mean(Xhat(:, :, :, X_test==1), 4));

[~, p] = ttest2(Xhat(:, :, :, X_test==0), Xhat(:, :, :, X_test==1), 'dim', 4);
p = squeeze(p);

pCorrect = squeeze(mean(class == repmat(permute(X_test, [1, 3, 4, 2]), [1, numT, numT, 1]), 4));

figure;
subplot(1, 3, 1); imagesc(time, time, m1-m0); eqClims; colorbar; axis xy; axis image;
subplot(1, 3, 2); imagesc(time, time, log10(p)); colorbar; axis xy; axis image;
subplot(1, 3, 3); imagesc(time, time, pCorrect); colorbar; axis xy; axis image;

%% 2.2.2. Decode using cross-validation
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

cfg = [];
cfg.folds = folds;
cfg.feedback = 'yes';
cfg.trainfun = 'train_array';
cfg.traincfg.feedback = 'yes';
cfg.traincfg.trainfun = 'train_LDA';
cfg.traincfg.traincfg.gamma = 0.05;
cfg.decodefun = 'decode_arrayGeneralization';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_LDA';

Xhat = decodeCrossValidation(cfg, X, Y);
class = (Xhat > 0);

% Plot
m0 = squeeze(mean(Xhat(:, :, :, X==0), 4));
m1 = squeeze(mean(Xhat(:, :, :, X==1), 4));

[~, p] = ttest2(Xhat(:, :, :, X==0), Xhat(:, :, :, X==1), 'dim', 4);
p = squeeze(p);

pCorrect = squeeze(mean(class == repmat(permute(X, [1 3 4 2]), [1, numT, numT, 1]), 4));

figure;
subplot(1, 3, 1); imagesc(time, time, m1-m0); eqClims; colorbar; axis xy; axis image;
subplot(1, 3, 2); imagesc(time, time, log10(p)); colorbar; axis xy; axis image;
subplot(1, 3, 3); imagesc(time, time, pCorrect); colorbar; axis xy; axis image;

%% 3. Decoding per sensor, using time points as featuers
%% 3.1. Using cross validation
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

cfg = [];
cfg.folds = folds;
cfg.feedback = 'yes';
cfg.trainfun = 'train_array';
cfg.traincfg.feedback = 'yes';
cfg.traincfg.trainfun = 'train_LDA';
cfg.traincfg.traincfg.gamma = 0.05;
cfg.decodefun = 'decode_array';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_LDA';

Xhat = decodeCrossValidation(cfg, X, permute(Y, [2 1 3]));
class = (Xhat > 0);

% Plot
m0 = squeeze(mean(Xhat(:, :, X==0), 3));
m1 = squeeze(mean(Xhat(:, :, X==1), 3));

[~, p] = ttest2(Xhat(:, :, X==0), Xhat(:, :, X==1), 'dim', 3);
p = squeeze(p);

pCorrect = squeeze(mean(class == repmat(permute(X, [1, 3, 2]), [1, length(label), 1]), 3));

figure;
subplot(1, 3, 1);
topoplot(label, m1 - m0); colorbar; eqClims;
subplot(1, 3, 2);
topoplot(label, log10(p)); colorbar;
subplot(1, 3, 3);
topoplot(label, pCorrect); colorbar;

figure; colormap(posneg(256));
topoplot(label, pCorrect); colorbar;
set(gca, 'clim', ((max(get(gca, 'clim')) - 0.5) * [-1 1]) + 0.5);

if (0)
    save('development/pCorrect_spat_2class.mat', 'pCorrect');
end


























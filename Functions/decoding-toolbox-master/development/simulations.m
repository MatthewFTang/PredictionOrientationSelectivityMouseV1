clear all;

% Add decoding functions to path
addpath('../');

%% 1. 2class problem
%% 1.1. Generate data
% Specificy parameters
numF = 2;
numN0 = 50;
numN1 = 150;

m0_true = [0 1]';          % True mean for class A
m1_true = [0 -1]';          % True mean for class B

S0_true = [1 .6; .6 1];               % True covariance for class A
S1_true = [1 .6; .6 1];               % True covariance for class B

% Generate data
Y0 = mvnrnd(repmat(m0_true', [numN0, 1]), S0_true)';
Y1 = mvnrnd(repmat(m1_true', [numN1, 1]), S1_true)';

Y = [Y0, Y1];

% Labels
X = [zeros(numN0, 1); ones(numN1, 1)];

%% 1.2. Performance analyses
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

% Do decoding using cross-validation
cfg = [];
cfg.trainfun = 'train_LDA';
cfg.traincfg.gamma = 0.05;
cfg.decodefun = 'decode_LDA';
cfg.folds = folds;
cfg.feedback = 'yes';
Xhat = decodeCrossValidation(cfg, X, Y)';

%% 1.3. Plot results
m0 = mean(Xhat(X==0));
m1 = mean(Xhat(X==1));

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
title(sprintf('mean diff: %.2f, m1log10(p): %.2f\npCorrect: %.2f', m1 - m0, log10(p), pCorrect));

% Biased?
[~, p] = ttest(Xhat)
pLeft = binocdf(sum(class==0), length(class), 0.5);
pRight = binocdf(sum(class==1), length(class), 0.5);
min([pLeft, pRight])*2










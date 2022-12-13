clear all;

addpath('../');

% Gives me: X (labels), Y (data), time, label
tmp = load('../data/testdata_orientation.mat');

Y = tmp.Y;              % Y: features x time x trials
X = tmp.X';
time = tmp.time;
label = tmp.label;

phi = X * (180/8);      % Presented orientation in degrees

numF = size(Y, 1);
numT = size(Y, 2);
numN = size(Y, 3);

%% Probabilistic classification
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, phi);

cfg = [];
cfg.folds = folds;
cfg.feedback = 'yes';
cfg.trainfun = 'train_array';
cfg.traincfg.feedback = 'yes';
cfg.traincfg.trainfun = 'train_probClass';
cfg.traincfg.traincfg.gamma = 0.01;
cfg.decodefun = 'decode_arrayGeneralization';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_probClass';

pPost = decodeCrossValidation(cfg, phi, Y);

% Extract classes
[~, class] = max(pPost, [], 1);
class = squeeze(class);

% Classification accuracy
correct = (class == repmat(permute((X+1), [1 3 2]), [numT, numT]));

[~, pCorrect] = ttest(correct*1, 1/8, 'dim', 3);
mCorrect = mean(correct, 3);

figure; colormap(jet(256));
subplot(2, 1, 1); imagesc(time, time, mCorrect); colorbar; eqClims(1/8); axis image; axis xy;
subplot(2, 1, 2); imagesc(time, time, log10(pCorrect) .* (pCorrect < 0.05)); colorbar; axis image; axis xy;

numC = 8;

pPostCorrect = zeros(numT, numT, numN);
for ic = 1:numC
    pPostCorrect(:, :, X==(ic-1)) = pPost(ic, :, :, X==(ic-1));
end

[~, ppPost] = ttest(pPostCorrect, 1/8, 'dim', 3);
mCorrect = mean(pPostCorrect, 3);

% Plot 
figure; colormap(jet(256));
subplot(2, 1, 1); imagesc(time, time, mCorrect); colorbar; eqClims(1/8); axis image; axis xy;
subplot(2, 1, 2); imagesc(time, time, log10(ppPost) .* (ppPost < 0.05)); colorbar; axis image; axis xy;

%% B&H model
% Create design matrix
numC = 8;

cfg = [];
cfg.numC = numC;
cfg.tuningCurve = 'vonMises';
cfg.kappa = 5;

design = designMatrix_BH(cfg, phi);

% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

cfg = [];
cfg.folds = folds;
cfg.feedback = 'yes';
cfg.trainfun = 'train_array';
cfg.traincfg.feedback = 'yes';
cfg.traincfg.trainfun = 'train_beamformer';
cfg.traincfg.traincfg.gamma = 0.01;
cfg.decodefun = 'decode_arrayGeneralization';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_beamformer';

Xhat = decodeCrossValidation(cfg, design, Y);

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = reshape(kernel*reshape(Xhat, [numC, numT*numT*numN]), [numT, numT, numN]);
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = exp(1i * (theta - repmat(permute(phi, [1 3 2]), [numT, numT, 1])) * (pi/180)*2);
r = abs(r) .* cos(angle(r));

[~, p] = ttest(r, 0, 'dim', 3);
mr = mean(r, 3);

% Plot 
figure; colormap(jet(256));
subplot(2, 1, 1); imagesc(time, time, mr); colorbar; eqClims; axis image; axis xy;
subplot(2, 1, 2); imagesc(time, time, log10(p) .* (p < 0.05)); colorbar; axis image; axis xy;














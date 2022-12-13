clear all;

% Add the toolbox
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

%% 1. No time dimension
sel_t = find(time >= .2, 1);

%% 1.1. Split into train- and test-set
cfg = [];
cfg.nFold = 2;
folds = createFolds(cfg, X);

X_train = X(folds{1});
Y_train = squeeze(Y(:, sel_t, folds{1}));
phi_train = phi(folds{1});

X_test = X(folds{2});
Y_test = squeeze(Y(:, sel_t, folds{2}));
phi_test = phi(folds{2});    

% Create design matrix
numC = 8;

cfg = [];
cfg.numC = numC;
cfg.tuningCurve = 'vonMises';
cfg.kappa = 5;

design = designMatrix_BH(cfg, phi_train);

% Train decoder
cfg = [];
cfg.gamma = 0.01;
decoder = train_pattern(cfg, design, Y_train);

% Decode
cfg = [];
cfg.gamma = 0.01;
Xhat = decode_pattern(cfg, decoder, Y_test);

% Shift accordign to orientation
Xhat_shifted = zeros(numC, length(folds{2}));

for ic = 1:numC
    Xhat_shifted(:, X_test==(ic-1)) = Xhat(mod((0:(numC-1)) + ic - 1 + floor(numC/2), numC) + 1, X_test==(ic-1));
end
    
% Average 
m_shifted = mean(Xhat_shifted, 2);
se_shifted = std(Xhat_shifted, [], 2) / sqrt(size(Xhat_shifted, 2));

m = zeros(numC, numC);
for ic = 1:numC
    m(:, ic) = mean(Xhat(:, X_test==(ic-1)), 2);
end

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = kernel*Xhat;
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - phi_test) * (pi/180)*2));
r = abs(r) * cos(angle(r));

mtheta = zeros(numC, 1);
for ic = 1:numC
    mtheta(ic) = mod(angle(mean(exp(1i * theta(X_test==(ic-1)) * (pi/180) * 2))), 2*pi) * (180/pi) * 0.5;
end

% Plot
figure;
subplot(2, 3, 1);
plot((0:(numC-1)) * (180/numC), m);
legend(arrayfun(@num2str, (0:(numC-1)) * (180/numC), 'uniformoutput', 0));
subplot(2, 3, 2); grid on;
boundedline((0:(numC-1)) * (180/numC) - 90, m_shifted, se_shifted);
subplot(2, 3, 3); hold on; grid on; axis equal;
plot(phi_test + 3*randn(1, length(folds{2})), theta, '.');
babline(0, 1, 'black');
plot((0:(numC-1)) * (180/numC), mtheta, '.red', 'markersize', 50);
title(sprintf('r: %.2f', r));

%% 1.2. Cross-validation
% Create folds
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

% Create design matrix
numC = 8;

cfg = [];
cfg.numC = numC;
cfg.tuningCurve = 'vonMises';
cfg.kappa = 5;

design = designMatrix_BH(cfg, phi);

% Decode using cross-validation
cfg = [];
cfg.feedback = 'yes';
cfg.folds = folds;
cfg.trainfun = 'train_pattern';
cfg.decodefun = 'decode_pattern';
cfg.decodecfg.gamma = 0.1;

Xhat = decodeCrossValidation(cfg, design, squeeze(Y(:, sel_t, :)));

% Shift accordign to orientation
Xhat_shifted = zeros(numC, numN);

for ic = 1:numC
    Xhat_shifted(:, X==(ic-1)) = Xhat(mod((0:(numC-1)) + ic - 1 + floor(numC/2), numC) + 1, X==(ic-1));
end
    
% Average 
m_shifted = mean(Xhat_shifted, 2);
se_shifted = std(Xhat_shifted, [], 2) / sqrt(size(Xhat_shifted, 2));

m = zeros(numC, numC);
for ic = 1:numC
    m(:, ic) = mean(Xhat(:, X==(ic-1)), 2);
end

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = kernel*Xhat;
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - phi) * (pi/180)*2));
r = abs(r) * cos(angle(r));

mtheta = zeros(numC, 1);
for ic = 1:numC
    mtheta(ic) = mod(angle(mean(exp(1i * theta(X==(ic-1)) * (pi/180) * 2))), 2*pi) * (180/pi) * 0.5;
end

% Plot
figure;
subplot(2, 3, 1);
plot((0:(numC-1)) * (180/numC), m);
legend(arrayfun(@num2str, (0:(numC-1)) * (180/numC), 'uniformoutput', 0));
subplot(2, 3, 2); grid on;
boundedline((0:(numC-1)) * (180/numC) - 90, m_shifted, se_shifted);
subplot(2, 3, 3); hold on; grid on; axis equal;
plot(phi + 3*randn(1, numN), theta, '.');
babline(0, 1, 'black');
plot((0:(numC-1)) * (180/numC), mtheta, '.red', 'markersize', 50);
title(sprintf('r: %.2f', r));

%% 2. Including temporal dimension
%% 2.1. Matching testing and decoding time
%% 2.1.1. Split into train- and test-set
cfg = [];
cfg.nFold = 2;
folds = createFolds(cfg, X);

X_train = X(folds{1});
Y_train = squeeze(Y(:, :, folds{1}));
phi_train = phi(folds{1});

X_test = X(folds{2});
Y_test = squeeze(Y(:, :, folds{2}));
phi_test = phi(folds{2});    

% Create design matrix
numC = 8;

cfg = [];
cfg.numC = numC;
cfg.tuningCurve = 'vonMises';
cfg.kappa = 5;

design = designMatrix_BH(cfg, phi);

% Train decoder
cfg = [];
cfg.feedback = 'yes';
cfg.trainfun = 'train_beamformer';
cfg.traincfg.gamma = 0.01;
decoder = train_array(cfg, design, Y_train);

% Decode
cfg = [];
cfg.feedback = 'yes';
cfg.decodefun = 'decode_beamformer';
Xhat = decode_array(cfg, decoder, Y_test);

% Shift accordign to orientation
Xhat_shifted = zeros(numC, numT, length(folds{2}));

for ic = 1:numC
    Xhat_shifted(:, :, X_test==(ic-1)) = Xhat(mod((0:(numC-1)) + ic - 1 + floor(numC/2), numC) + 1, :, X_test==(ic-1));
end
    
% Average 
m_shifted = mean(Xhat_shifted, 3);
[~, p_shifted] = ttest(Xhat_shifted, 0, 'dim', 3);

m = zeros(numC, numT, numC);
for ic = 1:numC
    m(:, :, ic) = mean(Xhat(:, :, X_test==(ic-1)), 3);
end

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = reshape(kernel*reshape(Xhat, [numC, numT*length(folds{2})]), [numT, length(folds{2})]);
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - repmat(phi_test, [numT, 1])) * (pi/180)*2), 2);
r = abs(r) .* cos(angle(r));

mtheta = zeros(numC, numT);
for ic = 1:numC
    mtheta(ic, :) = mod(angle(mean(exp(1i * theta(:, X_test==(ic-1)) * (pi/180) * 2), 2)), 2*pi) * (180/pi) * 0.5;
end

% Plot
figure;
for ic = 1:numC
    subplot(5, 3, ic);
    imagesc(time, 1:numC, m(:, :, ic)); colorbar;
end

subplot(5, 3, 10);
imagesc(time, 1:numC, m_shifted); colorbar;

subplot(5, 3, 13);
imagesc(time, 1:numC, log10(p_shifted)); colorbar;

subplot(5, 3, [11 12 14 15]); hold on; grid on; axis equal;
plot_t = find(time > .2, 1); 
plot(phi_test + 3*randn(1, length(folds{2})), theta(plot_t, :), '.');
babline(0, 1, 'black');
plot((0:(numC-1)) * (180/numC), mtheta(:, plot_t), '.red', 'markersize', 50);
title(sprintf('r: %.2f', r(plot_t)));

%% 2.1.2. Decode using cross-validation
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
cfg.traincfg.traincfg.gamma = 0.05;
cfg.decodefun = 'decode_array';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_beamformer';

Xhat = decodeCrossValidation(cfg, design, Y);

% Shift accordign to orientation
Xhat_shifted = zeros(numC, numT, numN);

for ic = 1:numC
    Xhat_shifted(:, :, X==(ic-1)) = Xhat(mod((0:(numC-1)) + ic - 1 + floor(numC/2), numC) + 1, :, X==(ic-1));
end
    
% Average 
m_shifted = mean(Xhat_shifted, 3);
[~, p_shifted] = ttest(Xhat_shifted, 0, 'dim', 3);

m = zeros(numC, numT, numC);
for ic = 1:numC
    m(:, :, ic) = mean(Xhat(:, :, X==(ic-1)), 3);
end

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = reshape(kernel*reshape(Xhat, [numC, numT*numN]), [numT, numN]);
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - repmat(phi, [numT, 1])) * (pi/180)*2), 2);
r = abs(r) .* cos(angle(r));

mtheta = zeros(numC, numT);
for ic = 1:numC
    mtheta(ic, :) = mod(angle(mean(exp(1i * theta(:, X==(ic-1)) * (pi/180) * 2), 2)), 2*pi) * (180/pi) * 0.5;
end

% Plot
figure;
for ic = 1:numC
    subplot(5, 3, ic);
    imagesc(time, 1:numC, m(:, :, ic)); colorbar;
end

subplot(5, 3, 10);
imagesc(time, 1:numC, m_shifted); colorbar;

subplot(5, 3, 13);
imagesc(time, 1:numC, log10(p_shifted)); colorbar;

subplot(5, 3, [11 14]); hold on; grid on; axis equal;
plot_t = find(time > .2, 1); 
plot(phi + 3*randn(1, numN), theta(plot_t, :), '.');
babline(0, 1, 'black');
plot((0:(numC-1)) * (180/numC), mtheta(:, plot_t), '.red', 'markersize', 50); axis equal;
title(sprintf('r: %.2f', r(plot_t)));

subplot(5, 3, 15); hold on; grid on;
plot(time, r);

if (0)
    r_time = r;
    save('development/r_time_orientation.mat', 'r_time');
end

%% 2.2. Temporal generalization
%% 2.2.1. Split into train- and test-set
cfg = [];
cfg.nFold = 2;
folds = createFolds(cfg, X);

X_train = X(folds{1});
Y_train = squeeze(Y(:, :, folds{1}));
phi_train = phi(folds{1});

X_test = X(folds{2});
Y_test = squeeze(Y(:, :, folds{2}));
phi_test = phi(folds{2});    

% Create design matrix
numC = 8;

cfg = [];
cfg.numC = numC;
cfg.tuningCurve = 'vonMises';
cfg.kappa = 5;

design = designMatrix_BH(cfg, phi);

% Train decoder
cfg = [];
cfg.feedback = 'yes';
cfg.trainfun = 'train_beamformer';
cfg.traincfg.gamma = 0.01;
decoder = train_array(cfg, design, Y_train);

% Decode
cfg = [];
cfg.feedback = 'yes';
cfg.decodefun = 'decode_beamformer';
Xhat = decode_arrayGeneralization(cfg, decoder, Y_test);

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = reshape(kernel*reshape(Xhat, [numC, numT*numT*length(folds{2})]), [numT, numT, length(folds{2})]);
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - repmat(permute(phi_test, [1 3 2]), [numT, numT, 1])) * (pi/180)*2), 3);
r = abs(r) .* cos(angle(r));

% Plot 
figure; colormap(posneg(256));
imagesc(time, time, r); colorbar; eqClims; axis image; axis xy;

%% 2.2.2. Decode using cross-validation
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

% Calculate covariance
numF = size(Y, 1);
numT = size(Y, 2);
numN = size(Y, 3);

S = zeros(numF, numF);
for it = 1:numT
    S = S + cov(squeeze(Y(:, it, :))')/numN;
end

gamma = 0.01;
S = (1-gamma)*S + gamma*eye(numF)*trace(S)/numF;
    
cfg = [];
cfg.folds = folds;
cfg.feedback = 'yes';
cfg.trainfun = 'train_array';
cfg.traincfg.feedback = 'yes';
cfg.traincfg.trainfun = 'train_pattern';
cfg.decodefun = 'decode_arrayGeneralization';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_pattern';
cfg.decodecfg.decodecfg.covariance = S;
%cfg.decodecfg.decodecfg.gamma = 0.1;

Xhat = decodeCrossValidation(cfg, design, Y);

% Extract single orientations and correlate
kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = reshape(kernel*reshape(Xhat, [numC, numT*numT*numN]), [numT, numT, numN]);
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - repmat(permute(phi, [1 3 2]), [numT, numT, 1])) * (pi/180)*2), 3);
r = abs(r) .* cos(angle(r));

% Plot 
figure; colormap(posneg(256));
imagesc(time, time, r); colorbar; eqClims; axis image; axis xy;

%% 3. Decoding per sensor, using time points as featuers
%% 3.1. Using cross validation
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
cfg.traincfg.traincfg.gamma = 0.05;
cfg.decodefun = 'decode_array';
cfg.decodecfg.feedback = 'yes';
cfg.decodecfg.decodefun = 'decode_beamformer';

Xhat = decodeCrossValidation(cfg, design, permute(Y, [2 1 3]));

% Extract single orientations and correlate
numS = length(label);

kernel = exp(1i * (0:(numC-1)) * (2*pi/numC));
Z = reshape(kernel*reshape(Xhat, [numC, numS*numN]), [numS, numN]);
theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation

r = mean(exp(1i * (theta - repmat(phi, [numS, 1])) * (pi/180)*2), 2);
r = abs(r) .* cos(angle(r));

% Plot 
figure; colormap(posneg(256));
topoplot(label, r); eqClims; colorbar;

if (0)
    r_spat = r;
    save('development/r_spat_orientation.mat', 'r_spat');
end




























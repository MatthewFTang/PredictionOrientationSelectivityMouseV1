clear all;

%% Load data
% Gives me: X (labels), Y (data), time, label
tmp = load('../data/testdata_orientation.mat');

Y = tmp.Y;              % Y: features x time x trials
X = tmp.X';
time = tmp.time;
label = tmp.label;

numF = size(Y, 1);
numN = size(Y, 3);

%% Single time opint
% Select time point
sel_t = find(time >= .2, 1);

% Split into train- and test-set
cfg = [];
cfg.nFold = 2;
folds = createFolds(cfg, X);

X_train = X(folds{1});
Y_train = squeeze(Y(:, sel_t, folds{1}));

X_test = X(folds{2});
Y_test = squeeze(Y(:, sel_t, folds{2}));

%% Train classifier
numC = length(unique(X));
gamma = 0.1;

% Calculate means and covariance matrix
m = zeros(numF, numC);
S = zeros(numF, numF);

for ic = 1:numC
    m(:, ic) = nanmean(Y_train(:, X_train==(ic-1)), 2);
    
    S = S + nancov(Y_train(:, X_train==(ic-1))');
end
S = S/numC;

% Regularize
S = (1-gamma)*S + gamma*eye(numF)*trace(S)/numF;

% Calculate weights and bias
W = S\m;
b = -0.5*diag(m'*W);

cfg = [];
cfg.gamma = 0.1;
cfg.demean = 'yes';
decoder = train_probClass(cfg, X_train, Y_train);

cfg = [];
cfg.demean = 'no';
pPost = decode_probClass(cfg, decoder, Y_test);

decoder = train_beamformer(cfg, X_train, Y_train);

%% Test classifier
p = exp(W'*Y_test + repmat(b, [1, length(X_test)]));
p = p ./ repmat(sum(p, 1), [numC, 1]);


mp = zeros(numC, numC);
for ic = 1:numC
    mp(:, ic) = mean(p(:, X_test==(ic-1)), 2);
end

figure;
imagesc(mp); colorbar; axis image;







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

%% 1. Select time point of interest
sel_t = find(time >= .12, 1);

%% Split into train- and test-set
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

%% Decode using standard LDA
gamma = 0.05;

S = 0.5*(cov(Y_train(:, X_train==0)') + cov(Y_train(:, X_train==1)'));
S = (1-gamma)*S + gamma*eye(numF)*trace(S)/numF;

m0 = mean(Y_train(:, X_train==0), 2);
m1 = mean(Y_train(:, X_train==1), 2);

w = 2 * (S\(m1 - m0))/((m1 - m0)'*(S\(m1 - m0)));

Xhat = w'*Y_test;
mean(Xhat(X_test==0))
mean(Xhat(X_test==1))

pCorrect = mean((Xhat > 0) == X_test)

w_ana = w;

%% Decode using gradient ascent LDA
cfg0 = [];

cfg0.log.stepI = 10;
cfg0.log.stepF = 5;

cfg0.numI = 10000;

cfg0.w0 = randn(numF, 1);

cfg0.alpha = 1;
cfg0.kappa = (10^(-4)/1)^(1/cfg0.numI);

cfg0.gamma = 0.05;
cfg0.lambda = .001;

cfg0.feedbackInt = 2;

decoder = train_LDA_iterative(cfg0, X_train, Y_train);

figure;
plot(decoder.log.w');

Xhat = decoder.w'*Y_test;
mean(Xhat(X_test==0))
mean(Xhat(X_test==1))

pCorrect = mean((Xhat > 0) == X_test)

w_GD = w;

%% Compare weights
figure;
plot(w_ana(:), w_GD(:), '.'); grid on; babline(0, 1);

figure; hold on;
plot(sort(abs(w_ana(:))), 'blue');
plot(sort(abs(w_GD(:))), 'green');
legend('ana', 'GD');


figure;
topoplot(label, w_GD);

%% Using cross-validation
cfg = [];
cfg.nFold = 5;
folds = createFolds(cfg, X);

% Decode using cross-validation
cfg = [];
cfg.feedback = 'yes';
cfg.folds = folds;
cfg.trainfun = 'train_LDA_iterative';
cfg.traincfg.log.stepI = 10;
cfg.traincfg.log.stepF = 5;
cfg.traincfg.numI = 10000;
cfg.traincfg.w0 = randn(numF, 1);
cfg.traincfg.alpha = 1;
cfg.traincfg.kappa = (10^(-4)/1)^(1/cfg.traincfg.numI);
cfg.traincfg.gamma = 0.05;
cfg.traincfg.lambda = .001;
cfg.traincfg.feedbackInt = 2;
cfg.decodefun = 'decode_LDA';

Xhat = decodeCrossValidation(cfg, X, squeeze(Y(:, sel_t, :)));

mean(Xhat(X==0))
mean(Xhat(X==1))

pCorrect = mean((Xhat > 0) == X)

%% 2. Select channel of interest
sel_chan = find(strcmp(label, 'MLO32'));

%% Split into train- and test-set
numN_train = 260;
numN_test = numN - numN_train;

i0 = find(X==0); i0 = i0(randperm(length(i0), numN_train/2));
i1 = find(X==1); i1 = i1(randperm(length(i1), numN_train/2));

i_train = [i0, i1];
i_test = setdiff(1:numN, i_train)';

X_train = X(i_train);
Y_train = squeeze(Y(sel_chan, :, i_train));

X_test = X(i_test);
Y_test = squeeze(Y(sel_chan, :, i_test));

%% Decode using standard LDA
gamma = 0.05;

S = 0.5*(cov(Y_train(:, X_train==0)') + cov(Y_train(:, X_train==1)'));
S = (1-gamma)*S + gamma*eye(numT)*trace(S)/numT;

m0 = mean(Y_train(:, X_train==0), 2);
m1 = mean(Y_train(:, X_train==1), 2);

w = 2 * (S\(m1 - m0))/((m1 - m0)'*(S\(m1 - m0)));

Xhat = w'*Y_test;
mean(Xhat(X_test==0))
mean(Xhat(X_test==1))

pCorrect = mean((Xhat > 0) == X_test)

w_ana = w;

%% Decode using gradient ascent LDA
cfg0 = [];

cfg0.log.stepI = 10;
cfg0.log.stepF = 5;

cfg0.numI = 10000;

cfg0.w0 = randn(numT, 1);

cfg0.alpha = 1;
cfg0.kappa = (10^(-4)/1)^(1/cfg0.numI);

cfg0.gamma = 0.05;
cfg0.lambda = .001;

cfg0.feedbackInt = 2;

decoder = train_LDA_iterative(cfg0, X_train, Y_train);

figure;
plot(decoder.log.w');

Xhat = decoder.w'*Y_test;
mean(Xhat(X_test==0))
mean(Xhat(X_test==1))

pCorrect = mean((Xhat > 0) == X_test)

w_GD = w;

%% Compare weights
figure;
plot(w_ana(:), w_GD(:), '.');

figure; hold on;
plot(sort(abs(w_ana(:))), 'blue'); grid on; babline(0, 1);
plot(sort(abs(w_GD(:))), 'green');
legend('ana', 'GD');

figure; hold on;
plot(time, w_ana, 'red');
plot(time, w_GD, 'blue');
legend('ana', 'GD');

%% 3. Both time and sensors
sel_time = 1:numT;%find(inrange(time, [0, .25]));
sel_chan = 1:numF;%find(ismember(label, ft_channelselection('M*O*', label)));

%% Split into train- and test-set
numN_train = 260;
numN_test = numN - numN_train;

i0 = find(X==0); i0 = i0(randperm(length(i0), numN_train/2));
i1 = find(X==1); i1 = i1(randperm(length(i1), numN_train/2));

i_train = [i0, i1];
i_test = setdiff(1:numN, i_train)';

X_train = X(i_train);
Y_train = Y(sel_chan, sel_time, i_train);
Y_train = reshape(Y_train, [length(sel_chan)*length(sel_time), numN_train]);

X_test = X(i_test);
Y_test = Y(sel_chan, sel_time, i_test);
Y_test = reshape(Y_test, [length(sel_chan)*length(sel_time), numN_test]);

%% Decode using standard LDA
%{
gamma = 0.05;

S = 0.5*(cov(Y_train(:, X_train==0)') + cov(Y_train(:, X_train==1)'));
S = (1-gamma)*S + gamma*eye(numT)*trace(S)/numT;

m0 = mean(Y_train(:, X_train==0), 2);
m1 = mean(Y_train(:, X_train==1), 2);

w = 2 * (S\(m1 - m0))/((m1 - m0)'*(S\(m1 - m0)));

Xhat = w'*Y_test;
mean(Xhat(X_test==0))
mean(Xhat(X_test==1))

pCorrect = mean((Xhat > 0) == X_test)

w_ana = w;
%}

%% Decode using gradient ascent LDA
numI = 10000;
alpha = 1;
kappa = (10^(-4)/1)^(1/numI);

gamma = 0.05;
lambda = .001;

A0 = Y_train(:, X_train==0);
A1 = Y_train(:, X_train==1);

m0 = mean(A0, 2);
m1 = mean(A1, 2);
m = 0.5*(m1 + m0);

numN0_train = sum(X_train==0);
numN1_train = sum(X_train==1);

A0 = A0 - repmat(m0, [1, numN0_train]);
A1 = A1 - repmat(m1, [1, numN1_train]);
A = Y_train - repmat(m, [1, numN_train]);

nu0 = sum(sum(A0.^2)) / (sum(X_train==0)-1) / (length(sel_chan)*length(sel_time));
nu1 = sum(sum(A1.^2)) / (sum(X_train==1)-1) / (length(sel_chan)*length(sel_time));
nu = sum(sum(A.^2)) / (numN_train-1) / (length(sel_chan)*length(sel_time));

w = randn((length(sel_chan)*length(sel_time)), 1);

log = [];
log.w = zeros((length(sel_chan)*length(sel_time)), numI);

tStart = tic;
for i = 1:numI
	Bw = 0.5*((1-gamma)*((w'*A0)*(A0'*w)/(numN0_train-1) + (w'*A1)*(A1'*w)/(numN1_train-1)) + gamma*(nu0 + nu1)*(w'*w));
    Bb = (1-gamma)*((w'*A)*(A'*w)/(numN_train-1)) + gamma*nu*(w'*w);    
   
    dw = (Bw.^2)\( ...
        ((1-gamma)/(numN_train)*A*(A'*w) + gamma*nu*w) ...
        *Bw - ...
        (0.5*((1-gamma)*(A0*(A0'*w)/(numN0_train-1) + A1*(A1'*w)/(numN1_train-1)) + gamma*(nu0 + nu1)*w)) ...
        *Bb) - lambda*w;
    
    w = w + alpha*dw;
       
    alpha = alpha*kappa;
    log.w(:, i) = w;

    if (toc(tStart) > 2)
        fprintf('Finished iteration %g/%g\n', i, numI);
        tStart = tic;
    end
end    

w = 2*w/((m1 - m0)'*w);

plotdata = log.w(randperm(length(sel_chan)*length(sel_time), 80), floor(linspace(0, numI-1, numI/100)) + 1);
%plotdata = log.w(randperm(length(sel_chan)*length(sel_time), 80), (end-1000):end);

figure;
plot(1:size(plotdata, 2), plotdata);

Xhat = w'*Y_test;
mean(Xhat(X_test==0))
mean(Xhat(X_test==1))

pCorrect = mean((Xhat > 0) == X_test)

%% Using train-function
sel_t = inrange(time, [.06, .15]);
sel_chan = ismember(label, ft_channelselection({'MZO*', 'MZF*'}, label));

% Split into train- and test-set
numN_train = 260;
numN_test = numN - numN_train;

i0 = find(X==0); i0 = i0(randperm(length(i0), numN_train/2));
i1 = find(X==1); i1 = i1(randperm(length(i1), numN_train/2));

i_train = [i0, i1];
i_test = setdiff(1:numN, i_train)';

X_train = X(i_train);
Y_train = squeeze(Y(sel_chan, sel_t, i_train));
Y_train = reshape(Y_train, [sum(sel_chan)*sum(sel_t), numN_train]);

X_test = X(i_test);
Y_test = squeeze(Y(sel_chan, sel_t, i_test));
Y_test = reshape(Y_test, [sum(sel_chan)*sum(sel_t), numN_test]);

w0 = randn(sum(sel_chan)*sum(sel_t), 1);

gamma = shrinkageGamma(Y_train, 1);

%% Train decoder
cfg = [];
cfg.log.stepI = 500;
cfg.log.stepF = 1;
cfg.numI = 10000;
cfg.w0 = w0;
cfg.alpha = .1;
cfg.kappa = 1;%(10^(-3)/1)^(1/cfg.numI);
cfg.gamma = gamma;
cfg.lambda_L1 = 0.01;
cfg.lambda_L2 = 0;
cfg.feedbackInt = 2;
decoder = train_LDA_iterative(cfg, X_train, Y_train);

%cfg = [];
%cfg.gamma = 0.05;
%decoder = train_LDA(cfg, X_train, Y_train);

[~, i] = sort(decoder.log.selF);
figure;
plot(1:decoder.log.nSelI, decoder.log.w(i, :));

Xhat = decode_LDA([], decoder, Y_test);

mean((Xhat > 0) == X_test)













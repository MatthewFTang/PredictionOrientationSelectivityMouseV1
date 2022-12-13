clear all;

addpath('~/-pimmos/decoding/');

% Gives me: X (labels), Y (data), time
tmp = load('../testdata_orientation.mat');

Y = tmp.Y;              % Y: features x time x trials
X = tmp.X';
time = tmp.time;
label = tmp.label;

numF = size(Y, 1);
numT = size(Y, 2);
numN = size(Y, 3);

phi = X * (180/8);      % Presented orientation in degrees

% Select time points, sensors, and collapse features
if (0)          % Manually
    selChan = find(ismember(label, ft_channelselection({'M*O*', 'M*P*'}, label)));
    selTime = find(inrange(time, [.06, .3]));

    fprintf('chan: %g, time: %g, memory: %.2f gb\n', length(selChan), length(selTime), (length(selChan)*length(selTime))^2*8/1024^3);
else            % Based on separate within-time and within-sensor decoding
    maxGB = 2;          % Maximum allowed size of covariance matrix, in GB
    
    numF = sqrt(maxGB*1024^3/8);
    
    load('r_spat_orientation.mat');
    load('r_time_orientation.mat');
    
    [~, i_spat] = sort(r_spat);
    [~, i_time] = sort(r_time);
    
    selChan = i_spat((end-round(sqrt(numF))+1):end);
    selTime = i_time((end-round(sqrt(numF))+1):end);
    
    fprintf('chan: %g, time: %g, memory: %.2f gb\n', length(selChan), length(selTime), (length(selChan)*length(selTime))^2*8/1024^3);
end    
    
Ycol = reshape(Y(selChan, selTime, :), [length(selChan)*length(selTime), numN]);

% Calculate gamma
tic
gamma = shrinkageGamma(Ycol, 1)
toc

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

% Decode using cross-validation
cfg = [];
cfg.feedback = 'yes';
cfg.folds = folds;
cfg.trainfun = 'train_beamformer';
cfg.traincfg.gamma = gamma;
cfg.decodefun = 'decode_beamformer';

Xhat = decodeCrossValidation(cfg, design, Ycol);

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
xlabel('True orientation + jitter');
ylabel('Decoded orientation');
title(sprintf('r: %.2f', r));


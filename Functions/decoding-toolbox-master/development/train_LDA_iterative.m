function decoder = train_LDA_iterative(cfg0, X, Y)

log = [];
log.cfg0 = cfg0;

A0 = Y(:, X==0);
A1 = Y(:, X==1);

m0 = mean(A0, 2);
m1 = mean(A1, 2);
m = 0.5*(m1 + m0);

numN = length(X);
numN0 = sum(X==0);
numN1 = sum(X==1);

numF = size(Y, 1);

A0 = A0 - repmat(m0, [1, numN0]);
A1 = A1 - repmat(m1, [1, numN1]);
A = Y - repmat(m, [1, numN]);

nu0 = sum(sum(A0.^2)) / (numN0-1) / numF;
nu1 = sum(sum(A1.^2)) / (numN1-1) / numF;
nu = sum(sum(A.^2)) / (numN-1) / numF;

alpha = cfg0.alpha;

w = cfg0.w0;

decoder = [];
decoder.log.nSelF = floor(numF/cfg0.log.stepF);
decoder.log.nSelI = floor(cfg0.numI/cfg0.log.stepI);
decoder.log.selF = randperm(numF, decoder.log.nSelF);
decoder.log.w = nan(decoder.log.nSelF, decoder.log.nSelI);

countLogI = 1;
iLog = 1;

tStart = tic;
for i = 1:cfg0.numI
	Bw = 0.5*((1-cfg0.gamma)*((w'*A0)*(A0'*w)/(numN0-1) + (w'*A1)*(A1'*w)/(numN1-1)) + cfg0.gamma*(nu0 + nu1)*(w'*w));
    Bb = (1-cfg0.gamma)*((w'*A)*(A'*w)/(numN-1)) + cfg0.gamma*nu*(w'*w);    
   
    dw = (Bw.^2)\( ...
        ((1-cfg0.gamma)/(numN-1)*A*(A'*w) + cfg0.gamma*nu*w) ...
        *Bw - ...
        (0.5*((1-cfg0.gamma)*(A0*(A0'*w)/(numN0-1) + A1*(A1'*w)/(numN1-1)) + cfg0.gamma*(nu0 + nu1)*w)) ...
        *Bb) - cfg0.lambda_L1*sign(w) - cfg0.lambda_L2*w;        

    w = w + alpha*dw;

    alpha = alpha*cfg0.kappa;
    
    if (countLogI == cfg0.log.stepI)
        decoder.log.w(:, iLog) = w(decoder.log.selF);
        iLog = iLog + 1;
        countLogI = 0;
    end
    countLogI = countLogI + 1;
       
    if (toc(tStart) > cfg0.feedbackInt)
        fprintf('Finished iteration %g/%g\n', i, cfg0.numI);
        tStart = tic;
    end
end    

decoder.W = (2*w/((m1 - m0)'*w))';
decoder.mY = 0.5*(m0 + m1)';

end

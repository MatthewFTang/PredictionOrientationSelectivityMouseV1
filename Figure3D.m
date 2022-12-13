clearvars
setUpPaths

%%
for len=1:23
    processData
    
    %%
    clear ExpectMeaned Sig
    
    trials=StoreType==1 & NewSupriseStore1==1;
    
    
    
    %%
    
    tempOris = StoreX2(trials);
    tempData = Store(trials,:,:);
    %% create the design matrix
    Orientations=0:30:150;
    X=tempOris;
    phi =Orientations(X);
    phi(phi==0)=180;
    
    % Generate design matrix
    fs=1/29.133;
    NewTime = -60*fs:fs:fs*90;
    
    TIO = NewTime >.25 & NewTime <1;
    Y = mean(tempData(:,:,TIO),3);
    
    %%
    modelChannels=8;
    n_ori_chans= modelChannels;
    
    funType = @(xx,mu) (cosd(xx-mu)).^(n_ori_chans-mod(n_ori_chans,2));
    
    xx = linspace(1,180,180);
    basis_set = nan(180,n_ori_chans);
    chan_center = linspace(180/n_ori_chans,180,n_ori_chans);
    
    for cc = 1:n_ori_chans
        basis_set(:,cc) = funType(xx,chan_center(cc));
    end
    
    stim_mask = zeros(length(phi),length(xx));
    
    for tt = 1:size(stim_mask,1)  % loop over trials
        stim_mask(tt,phi(tt))=1;
        
    end
    
    design = (stim_mask*basis_set)';
    numN=size(Y,1);
    
    
    Index =1:numN;
    cfg = [];
    cfg.nFold = 10;
    folds = createFolds(cfg, X);
    DataTime=squeeze(Y)';
    
    XhatTemp = nan(n_ori_chans,numN);
    for ii =1:length(folds)
        TestTrials =folds{ii};
        TrainTrials =find(~ismember(Index,TestTrials));
        Y_train=DataTime(:,TrainTrials);
        Y_test= DataTime(:,TestTrials);
        DesignTrain =design(:,TrainTrials);
        
        cfg = [];
        decoder=train_beamformerMT(cfg,DesignTrain,Y_train);
        XhatTemp(:,TestTrials) = decode_beamformer(cfg, decoder,Y_test);
    end
    Xhat=XhatTemp;
    
    
    
    
    
    %% re-aligns all trials on the presented orientation
    
    Xhat_algined = nan(size(Xhat),'single');
    uniPhi=unique(phi);
    numC=length(uniPhi);
    center=4;
    uniPhi=unique(phi);
    for ic = 1:length(uniPhi)
        ind = phi==uniPhi(ic);
        shiftInd= mod(ic-center:ic+center-1,8)+1;
        Xhat_algined(: ,ind)=...
            Xhat(shiftInd, ind);
    end
    plot(chan_center-110,mean(Xhat_algined,2))
    %%
    O =chan_center-110;
    decoded=zeros(1,size(Xhat,2));
    
    
    for trial =1:size(Xhat,2)
        d=Xhat_algined(:,trial);
        decoded(trial)=rad2deg( circ_mean(deg2rad(O*2)',d))/2;
        
    end
    edges=-90:10:90;
    [de bin]=histcounts(decoded,'BinEdges',edges);
    % decoded(perm,trial,:)=de;
    plot(bin(1:end-1),de)
    %%
    
    expectedOri = ExpectedOrientationsActual(trials);
    close all
    
    % supriseTrials = NewSupriseStore1(trials);
    
    presented = OrientationsActual(trials);
    error = presented-expectedOri;
    error = mod(error+90,180)-90;
    
    dat = Xhat_algined;
    datError =decoded;
    uniError = [-90 -60 -30 30 60];
    for ii=1:length(uniError)
        
        ind = error==uniError(ii);
        
        plot(mean(dat(:,ind),2),'LineWidth',2)
        hold on
        
        storeBig(len,ii,:)=mean(dat(:,ind),2);
        storeBigError(len,ii,:)=median(datError(:,ind),2);
        
    end
    legend(num2str(uniError'))
end


%%

lb = [-100 5 -10 -10];
ub = [100 180 10 10];
x0=[ 10 30 0 0];
x=-90:90;
close all
clear exp
f=@(c,xdata)c(3)+c(1).*(1/c(2)).*(xdata-c(4)).*exp(-(((xdata-c(4)).^2)/(2*c(2).^2)));



cols1 = linspecer(8,'green');
cols1(1:2,:)=[];

cols(1,:)=cols1(2,:);
cols(2,:)=cols1(4,:);
cols(3,:)=cols1(6,:);
cols(4,:)=cols1(3,:);

cols(5,:)=cols1(2,:);
cols(6,:)=cols1(1,:);
cols(3,:)=[];

close all
dat=squeeze(storeBigError(:,:,:));
datM=mean(dat);
datSE=std(dat)/sqrt(23);
[psy,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]  =  lsqcurvefit(f,x0,uniError',...
    datM',lb,ub);
fitLine = f(psy,x);
subplot(1,2,2)
plot(x,fitLine,'k','LineWidth',2)
hold on
for i=1:5
errorbar(uniError(i),datM(i),datSE(i),'o','CapSize',0,'MarkerFaceColor',cols(i,:),'Color',cols(i,:))
hold on
end
hold on
MakePrettyFigure;
xticks(-90:30:90)

subplot(1,2,1)
for i=1:5
    dat=squeeze(storeBig(:,i,:));
    datM=mean(dat);
    datSE=std(dat)/sqrt(23);
    
    errorbar(O,datM,datSE,'o','Color',cols(i,:),'CapSize',0,'MarkerFaceColor',cols(i,:),'DisplayName',num2str(uniError(i)))
    hold on
    
    
    y=double(datM);
    
    min_y=min(y);
    max_y=max(y);
    spread=(max_y-min_y);
    [max1, ind1] = max(y);
    
    lb = [-100 0 0 10];
    ub = double([1000 spread*1.15 180 100]);
    x0=double([min_y  spread O(ind1) 30]);
    
    psy1 = lsqcurvefit(@fitCic180,x0,O',y',lb,ub);
    fitLine1 = fitCic180(psy1,x);
    
    plot(x,fitLine1,'LineWidth',2,'Color',cols(i,:),'DisplayName',' ')
    
end
xticks(-90:90:90)
MakePrettyFigure
set(gcf,'POsition',[  277.6667  839.6667  822.6667  301.3333])
cd(figPath)
print(['DecodingBias' '.eps'],'-depsc')
subplot(1,2,1)
legend()

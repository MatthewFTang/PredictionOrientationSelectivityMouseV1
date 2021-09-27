clearvars
figPath ='/Volumes/SeagateMFT/QBI/Ehsan data/PredictionVisionExp/Ca Data/Figures'

dataPath=('/Volumes/SeagateMFT/QBI/Ehsan data/PredictionVisionExp/Ca Data/Suite2pResults/Saved')
dataPath=('E:\Saved\Saved')
figPath ='C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction\Figures'
cd(dataPath);
outputPath='C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction\Output'

addpath(genpath('/Users/u1093158/OneDrive - Australian National University/MATLAB/Toolboxes/vhlab-toolbox-matlab-master'))
fileNames = dir('*-SaveNoEpoch.mat');
%%
for len=1:length(fileNames)
    %%
    %%
cd(dataPath);
    
    load(fileNames(len).name)
    
    %%
    tic
    clear epochedFF
    TrialCount=1;
    epochedFF=nan(length(StoreType),sum(Iscell),161);
    for exp =1:length(StoreTimes)
        Ff=StoreFF{exp};
       % Ff=filtfast(Ff,1,[],'gaussian',2);
        TrialTimes = StoreTimes{exp};
        CaTimeStamp=StoreCaTimes{exp};
      for trial=1:length(TrialTimes)
                ind=find(CaTimeStamp'>TrialTimes(trial),1);
              
                temp =Ff(ind-60:ind+100,logical(Iscell))';
                
                epochedFF(TrialCount,:,:)=temp;
                
                TrialCount=TrialCount+1;
      end
            
    end
    toc
    
    %%
    Orientations=0:30:150;
    
    
    RandTrials = StoreType==2;
    Predict = StoreType ==1;
    ControlPredict = StoreType ==3;
    Expected = [1 diff(StoreX2(Predict))];
    Sequence = StoreX2(Predict);
    ExpectedA=abs(Expected);
    Store =epochedFF;
    
    AlwaysZeros = sum(sum(Store,1),3)==0;
         bad=any(any(isnan(Store),3),1);
                 Store= Store(:,~bad&~AlwaysZeros,:);

    Store=Store-mean(Store(:,:,30:60),3);
    TestOri=0:30:180;
    OrientationsStore=(StoreX2);
    OrientationsActual=Orientations(OrientationsStore);
    OrientationsActual(OrientationsActual==0)=180;
    FindConditionsPrediction
    
    %%
    if Nexps ==4
        ExpCond=3;
        
    else
        ExpCond=5;
    end
    
    %%
    clear ExpectMeaned Sig
    for expect =1:3
        if expect ==1
            trials =(RandTrials);
        elseif expect ==2
            trials=StoreType==1 & NewSupriseStore1==0;
        elseif expect ==3
            trials=StoreType==1 & NewSupriseStore1==1;
        elseif expect ==4
            trials=StoreType==3 & NewSupriseStore1==0;
        else
            trials=StoreType==3 & NewSupriseStore1==1;
        end
        
        for o=1:6
            BigStore{len,expect,o}=Store(trials&StoreX2==o,:,:);
        end
        
        %end
    end
end


%%
close all
clear S Ff Store StoreFF StoreCaTimes StoreSpikes epochedFF
for expect =[ 1 2 3]
    
   
    TC=1;
    Stacked=nan(7200,2000,161,'single');
    Ori =[];
    for o=1:6
        n=1;
        for m=1:23
            
            d= squeeze(BigStore{m,expect,o});
            Stacked(TC:TC+size(d,1)-1,n:n+size(d,2)-1,:)=d;
            n=n+size(d,2);
        end
        if expect ==1
            TC=TC+1200;
            Ori=[Ori ones(1,1200)*o];
            
        elseif expect ==2
            TC=TC+980;
            Ori=[Ori ones(1,980)*o];
            
        else
            TC=TC+180;
            Ori=[Ori ones(1,180)*o];
            
        end
    end
    Stacked=Stacked(1:TC-1,1:n-1,:);
    %%
    n=isnan(Stacked);
    n=~any(any(n,3),1);
    
    %% create the design matrix
    Orientations=0:30:150;
    X=Ori;
    phi =Orientations(X);
    phi(phi==0)=180;
    
    Y=Stacked(:,n,:);
    % Generate design matrix
    
    Y=permute(Y,[2 3 1]);
    Y=(Y(:,30:150,:));
    
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
    
    % Generate design matrix
    design = (stim_mask*basis_set)';
    numN=size(Y,3);
    numT=size(Y,2);
        for select =1:24  
            Xhat=nan(n_ori_chans,numT,numN,'single');

            select1 =randperm(size(Y,1))<=400;
            Index =1:numN;
            cfg = [];
            cfg.nFold = 10;
            folds = createFolds(cfg, X);
            Y1=Y(select1,:,:);
            %% does the IEM in one of two ways
            pim=1;
            if pim ==0
                for time = 1:numT
                    XhatTemp = nan(n_ori_chans,numN);
                    DataTime = squeeze(Y1(:,time,:));
                    for ii =1:length(folds)
                        
                        TestTrials =folds{ii};
                        TrainTrials =find(~ismember(Index,TestTrials))';
                        Y_train=DataTime(:,TrainTrials);
                        Y_test= DataTime(:,TestTrials);
                        DesignTrain =design(:,TrainTrials);
                        
                        W   = squeeze(Y_train*DesignTrain'*pinv(DesignTrain*DesignTrain')); %OLS
                        XhatTemp(:,TestTrials)=(pinv(W'*W)*W'*Y_test);
                    end
                    Xhat(:,time,:)  =  XhatTemp; %multiply test data by weights from training data
                    
                end
                
            else
                
                for time =1:numT
                    
                    DataTime = squeeze(Y1(:,time,:));
                    XhatTemp = nan(n_ori_chans,numN);
                    for ii =1:length(folds)
                        fprintf('Expect:%i | N:%i| select :%i | fold :%i\n',expect,time,select,ii)
                        TestTrials =folds{ii};
                        TrainTrials =find(~ismember(Index,TestTrials));
                        Y_train=DataTime(:,TrainTrials);
                        Y_test= DataTime(:,TestTrials);
                        DesignTrain =design(:,TrainTrials);
                        
                        cfg = [];
                        decoder=train_beamformerMT(cfg,DesignTrain,Y_train);
                        XhatTemp(:,TestTrials) = decode_beamformer(cfg, decoder,Y_test);
                    end
                    Xhat(:,time,:)=XhatTemp;
                end
            end
            
            %% multiples the IEM results by the basis set to go from channels to orientations
            reconstruction=nan(size(Xhat,3),size(Xhat,2),size(basis_set,1),'single');
            for time =1:size(Xhat,2)
                reconstruction(:,time,:)=squeeze(Xhat(:,time,:))'*basis_set';
            end
            reconstruction=permute(reconstruction,[3 2 1]);
            
            %
            %% re-aligns all trials on the presented orientation
            
            reconstruction_algined = nan(size(reconstruction),'single');
            uniPhi=unique(phi);
            numC=length(uniPhi);
            center=90;
            
            for ic = 1:180
                shiftInd= mod(ic-center:ic+center-1,180)+1;
                reconstruction_algined( :,:,phi==ic)=...
                    reconstruction(shiftInd,:, phi==ic);
            end
            
            S{expect,select}=reconstruction_algined;
        end
    
end
cd(outputPath)

save('IEMmodelAllNeuronsTime','-V7.3','S')

%%
cd(outputPath)
load('IEMmodelAllNeuronsTime.mat')

%%
close all
O=-90:89;
clear store
for e=[1:3]
    clear dat decoded
    for ii=1:24
        a=squeeze(S{e,ii});
        [e ii]
        parfor time =1:size(a,2)
            decoded=zeros(size(a,3),1);
        for n=1:size(a,3)
            d=a(:,time,n);
            decoded(n,:)=rad2deg( circ_mean(deg2rad(O*2)',d))/2;
            
        end
        edges=-90:10:90;
 [de bin]=histcounts(decoded,'BinEdges',edges);
        store(e,ii,time,:)=de;

        end
    end

end





%%
close all
fs=30;
NewTime = -60*fs:fs:fs*90;
NewTimeT=NewTime(30:150);
cols=linspecer(3);
x1=-90:90;
bins = -90:10:80;
bins=bins+5;
for e=1:3
TIO = NewTimeT>500 &NewTimeT<750;


    dat=squeeze(mean(store(e,:,TIO,:),3));
dat=dat./sum(dat,2);
        dM=squeeze(mean(dat))*100;
    dSE=squeeze(std(dat))/sqrt(24);
    
    errorbar(bins,dM,dSE,'o','Color',...
                cols(e,:),'MarkerFaceColor',cols(e,:),'CapSize',0)
    hold on
   
  %  plot([0 0],[.03 .09],'--k')
  
   spread =max(dM)-min(dM);
            [a ind1]=max(dM);
             lb = [-1 0 -90 10];
    ub = [20 20 180 100];
    x0=double([0  spread bins(ind1) 45]);
    
    
               psy =  lsqcurvefit(@fitCic180,x0,bins,...
            dM,lb,ub)
        fitLine=(fitCic180(psy,x1));
        plot(x1,fitLine,'Color',cols(e,:))
        
        
end
 box off
       set(gca,'FontSize',12,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom') 
    xticks(-90:90:90)
ylabel('Trials ( % )')
xlabel('Presented orientation\circ')
yticks(0:2:20)


set(gcf,'Position',[   680   505   278   473])
cd(figPath)
print('IEMmodelAllNeuronsTimeTuning.eps','-depsc')

%%
close all
clear psyStore
correct = abs(bins)<25;
for e=1:3
    for ii=1:size(store,2)
        for time =1:size(store,3)
    dM=squeeze((store(e,ii,time,:)))';
dM=dM/sum(dM);
dM=dM*100;
  psy=circ_r(deg2rad(bins'*2),dM');
        psyStore(e,ii,time,:)=psy;
        Decod(e,ii,time,:)=sum( dM(correct));
        end
    end
end
%%
var=2;
close all
clear hleng
for e=1:3
    
    
    %  histogram(decoded)
    dat=squeeze(psyStore(e,:,:));
    dM=squeeze(mean(dat));
    dSE=squeeze(std(dat))/sqrt(size(psyStore,1));
    
   hleng{e}= shadedErrorBar(NewTimeT,dM,dSE,{'Color',cols(e,:)},.5)
    hold on
    

end
xlim([-1000 2000])
xlabel('Time (ms)')
box off
set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')
set(gcf,'Position',[   469   516   351   470])
c=legend([hleng{1}.mainLine hleng{2}.mainLine hleng{3}.mainLine ],'Random','Expected','Unexpected');
set(gcf,'Renderer','painters')
cd(figPath)
print('IEMmodelAllNeuronsTime.eps','-depsc')

%%
%%
storeFilt=filtfast(store,3,[],'gaussian',1);
for e=1
    for ii=1:size(store,2)
        parfor time =1:size(store,3)
    dM=squeeze((storeFilt(e,ii,time,:)))';
dM=dM/sum(dM);
dM=dM*100;
   spread =max(dM)-min(dM);
            [a ind1]=max(dM);
             lb = [0 0 -5 5];
    ub = [20 50 5 180];
    x0=double([0  spread 0 30]);
    
    
               psy =  lsqcurvefit(@fitCic180,x0,bins,...
            dM,lb,ub);
        psyStore(e,ii,time,:)=psy;
        end
    end
end
        
%%




clearvars

dataPath=('E:\Saved\Saved')
figPath ='C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction\Figures'
cd(dataPath);
outputPath='C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction\Output'


fileNames = dir('*-SaveNoEpoch.mat');
%%
for len=1:23
    
    
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
    
    %% does the modelling
    n_Channels_Form=6; %number of channels in the form bank (range of pi radians) - must divide 180
    Sigma_Form=30;  %width of the channel (standard deviation of the Gaussian function)
    
    recoveryRatio=1.2;%1.2;
    close all
    a=figure;
    a.Color=[1 1 1];
    
    Oris=1:180;
    chanOris =rad2deg(pi/n_Channels_Form:pi/n_Channels_Form:pi);
    Channels_Form=zeros(n_Channels_Form,length(Oris));
    for ii=1:n_Channels_Form
        c=[0 1 chanOris(ii) Sigma_Form];
        
        Channels_Form(ii,:)=fitCic180(c,Oris);
    end
    
    OngoingAdapation=ones(1,size(Channels_Form,1));
    c=1;
    OngoingBasic1=ones(1,size(Channels_Form,1));
    
    clear testO adaptO Response_FormStore chanStore OngoingAdapationStore
    Channels_Form1=Channels_Form;
    chanStore=zeros(length(OrientationsActual),size(Channels_Form,1),size(Channels_Form,2));
    %ExpectedOrientationsActual(1)=180;
    cols = linspecer(6);
    expect=1;
    for expect_ratio=[0:.25:1];
        for o=30:30:180;
            c1=[0 expect_ratio o Sigma_Form];
            
            curve(o,:)= 1-fitCic180(c1,Oris);
        end
        adapt=1;
        for adapt_ratio = [0:.25:1];
                OngoingAdapation=ones(1,size(Channels_Form,1));

            for trial =1:length(ExpectedOrientationsActual)
                
                Test_Orientation=(OrientationsActual(trial));
                Expected_Orientation=(ExpectedOrientationsActual(trial));
                Channels_Form1=Channels_Form;
                
                
              %  Expected_Response=curve(Expected_Orientation,:);
             %   Channels_Form1=Channels_Form1.*Expected_Response;
                % end
                Channels_Form1=Channels_Form1.*OngoingAdapation';
                 Expected_Response=1-(Channels_Form1(:,Expected_Orientation))*expect_ratio;
                %
                              Channels_Form1=Channels_Form1.*Expected_Response;
 
                
                Response_Form=Channels_Form1(:,Test_Orientation);
                %  if StoreType(trial)==2
                OngoingAdapation= (1-(1-OngoingAdapation)/recoveryRatio);
                OngoingAdapation=OngoingAdapation.*abs(Response_Form.*adapt_ratio-1')';
                
                Response_FormStore(adapt,expect,trial,:)=Response_Form;
                
            end
            adapt=adapt+1;
        end
        expect=expect+1;
        
    end
    
    fs=29.133;
    NewTime = -60*fs:fs:fs*100;
    clear Bstore statStore
    close all
    %%
    size(Store)
    for expect =1:2
        
        if expect==1
            trials=StoreType==1;
        else
            trials=StoreType==2;
            
        end
        %% create the design matrix
        X=StoreX2(trials);
        phi =Orientations(X);
        phi(phi==0)=180;
        
        Y=Store(trials,:,:);
        Y1=Store(trials,:,:);
        
        Y=permute(Y,[2 3 1]);
        TIO = NewTime >250 & NewTime<800;
        Y=mean(Y(:,TIO,:),2);
        numN=size(Y,3);
        numT=size(Y,2);
        n_ori_chans= n_Channels_Form;
        Xhat=nan(n_ori_chans,numT,numN,'single');
        
        Index =1:numN;
        
        expectC=1;
        %%
        c= mod(randperm(sum(trials)),5)+1;
        
        Index= 1:sum(trials);
        
        %%
        
        clear r
        clear error Beta Rst CorrStore CorrStoreP
        for adapt=1:size(Response_FormStore,1)
            for expectC=1:size(Response_FormStore,2);
                design = squeeze(Response_FormStore(adapt,expectC,trials,:))';
                % design=cat(1,ones(size(design,2),1)',design);
                for  N=1:size(Y,1)
                    k = 5;
                    
                    
                     y =(squeeze(Y(N,:,:)));
                    X=design(:,:)';
                    
                    B= ridge(y,X,k,0);
                    %%
                   
                   % B=(pinv(X'*X))*X'*y;
                  % B=X\y
                       e=B(1) + B(2)*X(:,1) + B(3)*X(:,2) + B(4)*X(:,3) + B(5)*X(:,4)...
                      + B(6)*X(:,5) + B(7)*X(:,6);
                  %  e=X*(pinv(X'*X))*X'*y;
                    SSE = sum((e-y).^2);
                    SST = sum((mean(y)-y).^2);
                    
                   R2= 1-(SSE/SST);
                  % [B1,BINT,R,RINT,STATS] = regress(y,X);
                  % R
               
                    %%
                 %   x=0:.01:1;
                %    line = b(1) + b(2)*x + b(3)*x + b(4)*x +b(5)*x +b(6)*x +b(7)+x;
                %    plot(x,line)
                 %   plot(squeeze(Y(N,:,:)))
                    %%
                    
                    [RHO,PVAL] = corr(y,X);
                    
                    Beta(adapt,expectC,N,:)=B;
                                        Rst(adapt,expectC,N,:)=R2;
                                        CorrStore(adapt,expectC,N,:)=RHO;
                                        CorrStoreP(adapt,expectC,N,:)=PVAL;

                    
                end
            end
        end
        BigStore{expect,len}=Beta;
                BigStoreR{expect,len}=Rst;
                BigStoreCor{expect,len}=CorrStore;
                BigStoreCorP{expect,len}=CorrStoreP;

    end
end

%%
clear S SR  Scorr ScorrP
N=1;
for l=1:23
    for expect=1:2
        
        dat= BigStore{expect,l};
        Nnew=size(dat,3);
        S(:,:,expect,N:N+Nnew-1,:)=dat;
                SR(:,:,expect,N:N+Nnew-1)= BigStoreR{expect,l};
        Scorr(:,:,expect,N:N+Nnew-1,:)= BigStoreCor{expect,l};
        ScorrP(:,:,expect,N:N+Nnew-1,:)= BigStoreCorP{expect,l};

    end
    N=N+Nnew;
    
end

S=permute(S,[4 1 2 3 5]);
%SR=permute(SR,[4 1 2 3 5]);
Scorr=permute(Scorr,[4 1 2 3 5]);
ScorrP=permute(ScorrP,[4 1 2 3 5]);
%%

pref =squeeze(Scorr(:,1,1,1,:));
[ma maxI]=max(pref');
shiftCor=zeros(size(Scorr));
for ro=1:6
   ind = maxI==ro; 
    
   shiftCor(ind,:,:,:,:)=circshift(Scorr(ind,:,:,:,:),3-ro,5);
end
%%
cd(outputPath)
load('SigB.mat')
close all

for adapt=1:size(S,2)
    subplot(size(S,2),1,adapt)
    maxD=squeeze(mean(S(SigB,adapt,:,2,2:7),5));
    datM=mean(maxD);
    datSE=std(maxD)/sqrt(size(maxD,1));
    
    barwitherr(datSE,datM,'FaceColor',[.75 .75 .75])
    hold on
    for e=2:3
        [H,P,CI,STATS] =ttest(maxD(:,1),maxD(:,e))
        if P < .05
            plot(e,max(datM)*1.12,'*')
        end
    end
end
%%


n=1;
clear bestI
for N=1:size(S,1)
    d=Scorr(N,:,:,:,:);
    d1=squeeze(collapse(d,2:4));
    [i b]=max(d1);
    
    bestI(n,:,:,:)=squeeze(d(:,:,:,:,b));
    n=n+1
end
%%

uniOri =  deg2rad(unique(phi)*2);
for N=1:size(Scorr,1)
    for a=1:size(Scorr,2)
        for e=1:size(Scorr,2)
            for c=1:2
                dat= squeeze(S(N,a,e,c,2:7));
           %  plot(uniOri,dat)
          %   dat =[1 1 10000 1 1  1]'
             cirM(N,a,e,c,:)=   circ_r(uniOri',dat);
            end
        end
        
        
        
    end
end
    

%%
close all
for a=1:2
    dat =real(log(cirM(:,:,:,a)));
    subplot(2,1,a)
    dat=collapse(dat,1);
    imagesc(dat)
    colorbar
    
end
%%

for e=1:2
    subplot(2,1,e)
    d=S(SigB,:,:,e)*100;
    imagesc(collapse(d,1))
    colorbar
end
%%
    close all

for adapt =1:5
    subplot(5,1,adapt)
    dat=squeeze(Scorr(SigB,adapt,:,1,:));
    maxD=max(dat,[],3);
    datM=mean(maxD);
    datSE=std(maxD)/sqrt(size(dat,1));
    
    barwitherr(datSE,datM,'FaceColor',[.75 .75 .75])
    hold on
    for e=2:3
        [H,P,CI,STATS] =ttest(maxD(:,1),maxD(:,e))
        if P < .05
            plot(e,max(datM)*1.12,'*')
        end
    end
    
    
   % xticklabels({'Stimulus only','stimulus + model 0.25','stimulus + model 0.75'})
    xtickangle(45)
    set(gcf,'Position',[  446   593   254   389])
    box off
    ylabel('Beta')
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
end
cd(figPath)
print('BestBetaRudge','-depsc')

%%

cond=1;
  close all
adapt_ratios =0:.25:1;
for test =1:3
  ax(test)=  subplot(1,3,test )
  ST=S;
    if test ==1
    dat=squeeze(ST(SigB,:,1,cond,2:7));
    elseif test ==2
            dat=squeeze(ST(SigB,1,:,cond,2:7));

    else
            dat=squeeze(ST(SigB,4,:,cond,2:7));

    end 
    clear dat1
    for N=1:size(dat,1)
        
        dat1(N,:)=dat(N,:,maxI(SigB(N)));
        
    end

    datM=nanmean(dat1);
    datSE=nanstd(dat1)/sqrt(size(dat,1));
    
    barwitherr(datSE,datM,'FaceColor',[.75 .75 .75])
    set(gca,'XTickLabel',adapt_ratios)
    hold on
    
 
    %%
    %[p,tbl,stats] = ranova(dat1)
    xtickangle(45)
    set(gcf,'Position',[  446   593   254   389])
    box off
    if test ==1
        title('Adaptation Only')
            ylabel('Beta ')

    elseif test ==2
                title('Expectation Only')

    else
        
        title({['Expectation' ]
            ['Adaptation =0.25']})
        xlabel('Gain (a.u.)')
    end
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
    
end
match_ylim(ax)
set(gcf,'Position',[   446.0000  641.0000  687.5000  341.0000])
cd(figPath)
MakePrettyFigure
print('BestBetaRudge1','-depsc')
%%
cond=1;
cols = linspecer(8,'red')
cols =cols(3:2:end,:);
  close all
adapt_ratios =0:.25:1;
for test =1:3
  ST=S;
    if test ==1
    dat=squeeze(ST(SigB,:,1,cond,2:7));
    elseif test ==2
            dat=squeeze(ST(SigB,1,:,cond,2:7));

    else
            dat=squeeze(ST(SigB,4,:,cond,2:7));

    end 
    clear dat1
    for N=1:size(dat,1)
        
        dat1(N,:)=dat(N,:,maxI(SigB(N)));
        
    end

    datM=nanmean(dat1);
    datSE=nanstd(dat1)/sqrt(size(dat,1));
    
  h=errorbar(adapt_ratios,datM,datSE,'-o','Color',[cols(test,:)],...
        'MarkerFaceColor',[cols(test,:)],'CapSize',0,'LineWidth',2,...
        'MarkerSize',12)
    hold on


            ylabel('Beta ')

        xlabel('Gain (a.u.)')
    if test ==2
        plot([.1 .9],[.1, .1],'Color',cols(test,:),'LineWidth',3)
    elseif test ==3
                plot([.1 .9],[.105, .105],'Color',cols(test,:),'LineWidth',3)

    end
    
    
end
l=legend({'Adaptation Only','Expectation Only','','Adapation and Expectation'})
l.Box='off';
xlim([-.1 1.1])
set(l,'Position',[    0.1564    0.5916    0.3623    0.1369])
set(gcf,'Position',[   446.0000  491.5000  863.5000  414.5000
])
yticks(0.02:.04:1)
cd(figPath)
MakePrettyFigure
print('BestBetaRudge1','-depsc')
%%
dat2=[];
for test =1:3
  ax(test)=  subplot(1,3,test )
  ST=S;
    if test ==1
    dat=squeeze(ST(SigB,:,1,cond,2:7));
    elseif test ==2
            dat=squeeze(ST(SigB,1,:,cond,2:7));

    else
            dat=squeeze(ST(SigB,4,:,cond,2:7));

    end 
    clear dat1
      for N=1:size(dat,1)
        
        dat1(N,:)=dat(N,:,maxI(SigB(N)));
        
      end
    
    dat2=[dat2 dat1]
end

   
    dat  =[[1:15] ;dat2]
        csvwrite('OverallFits.csv',dat)

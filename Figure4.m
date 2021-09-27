clearvars
dataPath=('E:\Saved\Saved\Anaesthesia')
figPath ='C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction\Figures'
cd(dataPath);
outputPath='C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction\Output'


%%
fileNames = dir('*-SaveEpochedNew.mat');

for len=1:length(fileNames)
    %%
cd(dataPath)
    
    load(fileNames(len).name)
    
    
    %%
    
    RandTrials = StoreType==2;
    Predict = StoreType ==1;
    ControlPredict = StoreType ==3;
    Expected = [1 diff(StoreX2(Predict))];
    Sequence = StoreX2(Predict);
    ExpectedA=abs(Expected);
    
    AlwaysZeros = sum(sum(Store,1),3)==0;
    
    Store= StoreRaw(:,~AlwaysZeros,:);
    Store=Store-mean(Store(:,:,30:60),3);
    TestOri=0:30:180;
    OrientationsStore=(StoreX2);
   FindConditionsPrediction
    %%
    
    
    %%
    if Nexps ==4
        ExpCond=3;
        
    else
        ExpCond=5;
    end
    %%
    fs=29.133;
    
    % cells = find(Iscell);
    clear ExpectMeaned Sig fano ExpectSE
    for expect =1:ExpCond
        
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
        tempOris = StoreX2(trials);
        length(tempOris)
        
        tempStore = Store(trials,:,:);
        clear Expectataion
        for ori =1:6
            d=(tempStore(tempOris==ori,:,:));
            %    sum(tempOris==ori)
            Expectataion(1:size(d,1),:,:,ori)=d-mean(d(:,:,1:60),3);
            ExpectMeaned(:,:,ori,expect)=collapse(d,1);
                        ExpectSE(:,:,ori,expect)=squeeze(std(d,1)/sqrt(size(d,1)));

        end
        
        %%
        NewTime = -1/fs*60:1/fs:1/fs*100;
        
        close all
        
        %%
        
        Act= squeeze(tempStore(:,:,:));
        ActBase= mean(Act(:,:,30:60),3);
        ActEvoked= mean(Act(:,:,60:90),3);
        clear actualE PermE
        Nperms=10000;
        PermE=zeros(Nperms,6,size(Act,2));
        for o=1:6
            actualE(o,:)=   mean(ActEvoked(tempOris==o,:),1)-mean(ActBase(tempOris==o,:),1);
            
            parfor perm=1:Nperms
                if mod(perm,100)==0
                    fprintf('Date:%2.0f | expect:%2.0f |Ori:%2.0d|perm:%2.0d \n',len,expect,o,perm)
                end
                order = randperm(length(tempOris));
                Pori=tempOris(order);
                PermE(perm,o,:)=   mean(ActEvoked(Pori==o,:))-mean(ActBase(Pori==o,:));
                
                
            end
        end
        %%
        for N=1:size(actualE,2)
            for o=1:6
                
                if actualE(o,N)<0
                    a=sum(PermE(:,o,N)<actualE(o,N));
                    
                else
                    a=sum(PermE(:,o,N)>actualE(o,N));
                end
                Sig(expect,N,o)=a/Nperms;
                
            end
            
        end
        
        SigNeuorns =Sig<.025;
        %%
        
        for o=1:6
            d= mean(Act(tempOris==o,:,60:75),3);
            dVar = var(d);
            dMean=mean(d);
            fano(:,expect,o)=dVar./dMean;
        end
        %%
        Ori =0:30:150;
        close all
        s=find(any(squeeze(SigNeuorns(expect,:,:))'));
    end
    cd(outputPath)
    save([fileNames(len).name(1:end-19) '-TuningSignifcanceTestingNewSortingAna'],'ExpectMeaned','Sig','fano','ExpectSE')
    
end
%%

    cd(outputPath)
td=dir('*TuningSignifcanceTestingNewSortingAna.mat');


OverallExpectedC=[];
SigNeuronsOverallC=[];
OverallExpectedC=nan(1000,161,6,6);
OverallExpectedSE=nan(1000,161,6,6);

SigNeuronsOverallC=nan(6,1000,6);
FanoStore=[];
c=1;
for i=1:length(td)
    load(td(i).name)
    %  T=OverallExpected{i};
    %  loa
    td(i).name
    N=size(ExpectMeaned,1);
    OverallExpectedC(c:c+N-1,:,:,1:size(ExpectMeaned,4))=ExpectMeaned;
        OverallExpectedSE(c:c+N-1,:,:,1:size(ExpectMeaned,4))=ExpectSE;

    SigNeuronsOverallC(1:size(Sig,1),c:c+N-1,:)=Sig;
    Control(c:c+N-1)=size(ExpectMeaned,4);
    FanoStore(c:c+N-1,1:size(ExpectMeaned,4),:)=fano;
    SessStore(c:c+N-1,:)=i;
    c=c+N;
end
%%
    cd(outputPath)
load('SigPredictionCirc.mat');


OverallExpectedAwake=[];
SigNeuronsOverallAwake=nan(6,1000);
c=1;
for i=1:size(StoreExpect,2)
    ExpectMeaned =StoreExpect{i};
    N=size(ExpectMeaned,1);
    datP=StoreStats{i};
    OverallExpectedAwake(c:c+N-1,:,:,1:size(ExpectMeaned,4))=ExpectMeaned;
    SigNeuronsOverallAwake(1:size(datP,1),c:c+N-1,:)=datP;
    c=c+N;
end
c
%%
Ncount=1;
clear fitStore diffStore
TIO = NewTime>0 & NewTime<.6;
Curve= collapse(OverallExpectedC(:,TIO,:,:),2);

for N=Sig
    if Control(N)==3
        t=[1 2 3];
    else
        t=[1:5];
        
    end
    for expect =[ t  ]
        dat = squeeze(Curve(N,:,expect));
        
        y=dat;
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        
        lb = [-100 -1 0 15];
        ub = [10000 spread*1.15 180 100];
        x0=double([min_y  spread Ori(ind1) 30]);
        
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub,options);

        fitStore(Ncount,expect,:)=psy1;
        
    end
    ControlCounter(Ncount)=expect;
    Ncount=Ncount+1;
end
%%
%%
Ncount=1;
clear fitStoreAwake 
TIO = NewTime>0 & NewTime<.6;
CurveAwake= collapse(OverallExpectedAwake(:,TIO,:,:),2);
SigNeuornsAwake=SigNeuronsOverallAwake<.05;

S2=(squeeze(SigNeuornsAwake([1 3],:,:)));
SigAwake=find((any(S2)))
for N=SigAwake
    t=1:3;
    for expect =[ t  ]
        dat = squeeze(CurveAwake(N,:,expect));
        
        y=dat;
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        
        lb = [-100 -1 0 15];
        ub = [10000 spread*1.15 180 100];
        x0=double([min_y  spread Ori(ind1) 30]);
        
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub,options);

        fitStoreAwake(Ncount,expect,:)=psy1;
        
    end
    ControlCounter(Ncount)=expect;
    Ncount=Ncount+1;
end
Ncount
%%

var=2;
close all
hold on

good = abs(fitStore(:,1,3)-fitStore(:,1,3))<180;
sum(good)
for compare =1:3
    if compare ==1
        dat= fitStore(:,[1 3],var);
    elseif compare ==2
                dat= fitStore(:,[2 3],var);

    else
                     dat= fitStore(:,[1 2],var);
   
    end
    dat=log(dat);
   ax(compare)= subplot(1,3,compare);
        plot(dat(:,1),dat(:,2),'k.',...
            'MarkerFaceColor','k','MarkerSize',12)
        hold on
        plot([-2 7],[-2 7],'k')

axis square

[H,P,CI,stats] =ttest(dat(:,1),dat(:,2));
title(sprintf('N = %2.0f | t = %2.2f, p = %2.4f',...
    size(fitStore,1),real(stats.tstat),real(P)))
box off
real(stats.tstat),
if compare ==3
    ylabel('Amplitude expected')

else
ylabel('Amplitude unexpected')
end
if compare ==2
    xlabel('Amplitude expected')

else
xlabel('Amplitude random')

end
end

match_ylim(ax);
match_xlim(ax)
% end
set(gcf,'Position',[   103   443   930   663])
cd(figPath)
print(gcf, 'ScatterPrectvsRandAna.pdf')



%%
close all
var=2;
    Ana=log(fitStore(:,[1  3],var));
AwaSuprise =Ana(:,2)-Ana(:,1);


    Awake=log(fitStoreAwake(:,[1  3],var));
AwakeSuprise =Awake(:,2)-Awake(:,1);
datM= [mean(AwaSuprise) mean(AwakeSuprise)];
datSE=[std(AwaSuprise)/sqrt(length(AwaSuprise)) std(AwakeSuprise)/sqrt(length(AwakeSuprise))];
    [H,P,CI,stats] =ttest2(real(AwaSuprise),real(AwakeSuprise))

 barwitherr(datSE,datM,'FaceColor',[.7 .7 .7],'basevalue',Bv)
    set(gca,'XTickLabel',{'Anaesthetized','Awake'})
    xlim([.5 2.5])
    ylabel('Suprise size (Expected minus Random)')
yticks(0:.2:1)
        xtickangle(45)
    box off
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom','FontName','Helvetica')
    set(gcf,'Position',[   680   653   261   445])
title(sprintf('N = %2.0f | t = %2.2f, p = %2.4f',...
    size(fitStore,1),real(stats.tstat),real(P)))
cd(figPath)
MakePrettyFigure
print(gcf, 'AwakeVsAna.pdf')

%%
Ncount=1;
clear Shifted
for N=Sig
    for expect =[ 1:5]
        dat = squeeze(Curve(N,:,expect));
        
        y=dat;
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        t=squeeze(OverallExpectedC(N,:,:,expect));
        Shifted(Ncount,expect,:,:)=circshift(t,4-ind1,2);
    end
    
    
    
    Ncount=Ncount+1;
end
%%
close all
Oris=(0:30:150)-90;

cols = linspecer(6);
for e=1:3
    subplot(3,1,e)
    for o=1:6
        dat= squeeze(Shifted(:,e,:,o));
        datM= squeeze(nanmean(dat));
        datSE = squeeze(nanstd(dat))/sqrt(length(dat));
        h(o)=shadedErrorBar(NewTime,datM,datSE,{'Color',cols(o,:)},.5);
        %plot(NewTime,datM)
        xlim([-.7 2])
        ylim([-10 10])
        yticks(-10:10:20)
        
        hold on
        
        
    end
    title(size(Shifted,1))
    
    plot([-20 20],[0 0],'k')
    plot([0 0],[-100 100],'k')
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
    box off
    ylim([-4 10])
end
legend(num2str(Oris'))
c=legend([h.mainLine  ],num2str(Oris'));
set(gcf,'Position',[680   564   216   534])
set(gcf,'Renderer','painters')
cd(figPath)
MakePrettyFigure
print(['AlignedNeuronsAna' '.eps'],'-depsc')


%%
close all
X1=-90:90
cols = linspecer(6);
for e=1:3
    dat= squeeze(mean(Shifted(:,e,TIO,:),3));
    datM= squeeze(mean(dat));
    datSE = squeeze(std(dat))/sqrt(length(dat));
    errorbar(Oris,datM,datSE,'o','MarkerFaceColor',cols(e,:),...
        'MarkerSize',6,'CapSize',0,'MarkerEdgeColor',cols(e,:))
    % xlim([-.7 2])
    xlim([-90 90])
    %  ylim([-8 15])
    %yticks(-10:10:20)
    
    hold on
    
    
    y=datM;
    min_y=min(y);
    max_y=max(y);
    spread=(max_y-min_y);
    [max1, ind1] = max(y);
    
    lb = [-100 min_y*1.2 0 -100];
    ub = [100 spread*1.15 180 100];
    x0=double([min_y  spread Ori(ind1)-90 20]);
    
    
    
    psy1 = lsqcurvefit(@fitCic180,x0,Oris,datM,lb,ub,options)
    xticks(0:90:180)
    fitLine1 = fitCic180(psy1,X1);
    
    plot((X1),fitLine1,'-','Color',cols(e,:))
    
    
    
    
    title(size(Shifted,1))
    
    %   plot([-20 20],[0 0],'k')
    %     plot([0 0],[-100 100],'k')
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
    
    
end
set(gcf,'Position',[680   564   216   216])

box off
legend('Random','','Expected','','Unexpected')
cd(figPath)
MakePrettyFigure
print(['AlignedNeuronsCurve' '.eps'],'-depsc')


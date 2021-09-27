clearvars
overallPath = 'C:\Users\Desktop\OneDrive - Australian National University\Current work\Rodent V1 and prediction'
dataPath=('E:\Saved\Saved')
outputPath=[overallPath '\Output']
figPath =[overallPath '\Figures']
cd(dataPath);
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
    
    RandTrials = StoreType==2;
    Predict = StoreType ==1;
    ControlPredict = StoreType ==3;
    Expected = [1 diff(StoreX2(Predict))];
    Sequence = StoreX2(Predict);
    ExpectedA=abs(Expected);
    Store= epochedFF(:,:,:);
    
    AlwaysZeros = sum(sum(Store,1),3)==0;
    Store= epochedFF(:,~AlwaysZeros,:);
    Store=Store-mean(Store(:,:,30:60),3);
    TestOri=0:30:180;
    OrientationsStore=(StoreX2);
    FindConditionsPrediction
   
    %%
    if Nexps ==4
        ExpCond=3;
    else
        ExpCond=5;
    end
    %%
    fs=29.133;
    clear ExpectMeaned Sig fano pStore sCir pStore pStoreTime  Estore pStoreTimeError Expectataion
    for expect =1:ExpCond
        
        NewTime = -1/fs*60:1/fs:1/fs*100;

        TIO = NewTime >0.250 & NewTime <1.000;
    
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
        for ori =1:6
            d=(tempStore(tempOris==ori,:,:));
            size(d)
            %    sum(tempOris==ori)sCir
            d= d-mean(d(:,:,1:60),3);
            d1=mean(d(:,:,TIO),3);
            Expectataion(:,ori,expect)=std(d1)/sqrt(size(d,1));
            ExpectMeaned(:,:,ori,expect)=collapse(d,1);
            
        end
        
        
        
        tempStoreNoTime = mean(tempStore(:,:,TIO),3);
        for N=1:size(tempStore,2)
            [p h comp]= myanova(tempStoreNoTime(:,N),tempOris);
            pStore(expect,N,:)=p;
        end
        
        for N=1:size(tempStore,2)
            for time=1:161
                pStoreTime(expect,N,time,:)=myanova(tempStore(:,N,time),tempOris);
            end
        end
        %%
        O=0:30:150;
        for N=1:size(ExpectMeaned,1)
            dat =squeeze(ExpectMeaned(N,:,:,expect));
            dat=dat-min(min(dat));
            for t=1:161
                
                r =circ_r(deg2rad(O*2)',dat(t,:)');
                sCir(expect,N,t)=r;
            end
        end

    end
    %%
    clear St sigMod
    if ExpCond==3
        St(RandTrials)=1;
        
        St(StoreType==1 & NewSupriseStore1==0)=2;
        
        St(StoreType==1 & NewSupriseStore1==1)=3;
    else
        
        St(RandTrials)=1;
        
        St(StoreType==1 & NewSupriseStore1==0)=2;
        
        St(StoreType==1 & NewSupriseStore1==1)=3;
        St(StoreType==3 & NewSupriseStore1==0)=4;
        St(StoreType==3 & NewSupriseStore1==1)=5;
        
        
    end
    
    
    StoreNoTime=collapse(Store(:,:,TIO),3);
    for N=1:size(Store,2)
        ExpectMeanedNo=squeeze(collapse(ExpectMeaned(N,TIO,:,[1 3]),[2 4]));
        [best bestI]=max(ExpectMeanedNo)
        ind = StoreX2==bestI;
        ind2=St==1|St==3;
        storeTypeInd = St(ind&ind2);
        
        dat= StoreNoTime(ind&ind2,N);
        
        [p h comp]= myanova(dat,storeTypeInd');
        sigMod(N,:)=p;
    end
    
    %%
    StoreCirc{len}=sCir;
    StoreStats{len}=pStore;
    StoreExpect{len}=ExpectMeaned;
    StoreMod{len}=sigMod;
    StoreStatTime{len}=pStoreTime;
    StoreExpectation{len}=Expectataion;
    
end
cd(outputPath)
save('SigPredictionCirc','StoreMod','StoreExpect','StoreStats')
%%
nCount=1;
BigS=[]
BigStat=[];
BigExpect=[]'
BigStatTime=[];
BigMod=[];
BigExpectation=[];
Control=[];
N=1;
for len=1:length(StoreCirc)
    dat =StoreCirc{len};
    datP=StoreStats{len};
    newN=size(dat,2);
    
    BigS(1:size(dat,1),N:N+newN-1,:)=dat;
    BigStat(1:size(dat,1),N:N+newN-1)=datP;
    BigStatTime(1:size(dat,1),N:N+newN-1,:)=StoreStatTime{len};
    
    BigExpect(N:N+newN-1,:,:,1:size(dat,1))=StoreExpect{len};
    BigExpectation(N:N+newN-1,:,1:size(dat,1))=StoreExpectation{len};
    
    BigMod(N:N+newN-1)=StoreMod{len};
    Control(N:N+newN-1)=size(dat,1);
    Session(len)=size(dat,1)
    N=N+newN;
end

%%
SigB=find(BigMod<.05);
cd(outputPath)

save('SigB','SigB')
%%
close all
sigN=find(any(BigStat([1 3],:)<.05));
sprintf('N :%i orientation selective neurons',length(sigN))
TIO = NewTime >.25 & NewTime <1;

dat= mean(BigS([1 3],:,TIO),3);
B=dat(2,:)-dat(1,:);
[q bB]=sort(B);
histogram(B)
sprintf('N :%i response to preferred orientation modulated by expectation',length(SigB))

%%
clear T T1
close all
nc=1;
for N=sigN
    
    dat =squeeze(mean(BigExpect(N,TIO,:,1:3),2));
    T(nc,:)=max(dat);
    T1(nc,:)=BigMod(N);
    
    nc=nc+1;
end
T=log(T*100);
hold on
axis square
box off
plot([-1 4],[-1 4],'k')
for i=1:length(T)
    if T1(i)<.05
        col='r';
    else
        col='k';
    end
    
    plot((T(i,1)),(T(i,3)),'o','MarkerFaceColor',col,...
        'Color',col,'MarkerSize',5)
    
end
ylabel('Response Unexpected (log dF/F)')
xlabel('Response Random (log dF/F)')
set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')
title(sprintf('Sig neurons :%2.0f \ %2.0f',sum(T1<.05),length(T)))
cd(figPath)
print(['Response Preferred Modulated Neurons' '.eps'],'-depsc')
%%
close all
ResponseSigMod =[];
c=1;
for N=1:size(BigExpect,1)
    ExpectMeanedNo=squeeze(collapse(BigExpect(N,TIO,:,[1 3]),[2 4]));
    [best bestI]=max(ExpectMeanedNo)
    dat =squeeze(collapse(BigExpect(N,TIO,bestI,[1 3]),[2 ]));
    ResponseSigMod(c,:)=dat;
    c=c+1;
end

%%
close all

dif= (ResponseSigMod(:,2)-ResponseSigMod(:,1))*100;
histogram(dif,30,'FaceColor',[.75 .75 .75])
sum(dif>0)
ylabel('Neurons')
xlabel('Reponse Unexpected minus Random (dF/F %)')
set(gca,'FontSize',12,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')
axis square
box off
title(['N= ' num2str(length(dif))])
set(gcf,'Renderer','painters')
cd(figPath)
print(['ModulationByPredictionN' '.eps'],'-depsc')

%%
close all
b=-1:.25:5;
cols = linspecer(3)
histogram(T(:,1),b,'FaceColor',cols(1,:),'Normalization','probability')
hold on
histogram(T(:,3),b,'FaceColor',cols(3,:),'Normalization','probability')
xlabel('log dF/F')
ylabel('Neurons (P)')
axis square
box off
legend('Random','Unexpected')
set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')
cd(figPath)
print(['Response Preferred Modulated Neurons Histogram' '.eps'],'-depsc')




%%
close all
[sortR sortI]=sort(BigStat(3,:)');
N=sortI(6)
X1=0:180;
Ori=0:30:150;
TIO = NewTime >.25 & NewTime <1;

dat1=squeeze(BigExpect(N,:,:,1:3))*100;

mMin=min(min(min(dat1)));
mMax=max(max(max(dat1)));
dif =mMax-mMin;
%legend('R','expected','UnExpceted')
count=1;
clear ax
for i=1:6
    ax(i)= subplot(1,7,i);
    
    fh=fill([0.25 ,0.25 ,1,1],[mMax mMin mMin mMax],[.85 .85 .85],'EdgeColor','none');
    
    for e=[1:3]
        dat=squeeze(BigExpect(N,:,:,e))*100;
        dat = filtfast(dat,1,[],'gaussian',1);
        % plot([0 0],[mMin mMax],'k','LineWidth',1)
        
        hold on
        
        plot(NewTime,dat(:,i),'Color',cols(e,:),'LineWidth',2)
        
        plot([NewTime(1) NewTime(end)],[0 0 ],'k','LineWidth',1)
        plot([0 0],[mMin mMax],'k','LineWidth',1)
        
        
        xlim([-1 2.0])
        
        
        axis off
        box off
        set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
            'YAxisLocation','left','LineWidth',1.5,...
            'XMinorTick','off','YMinorTick','off'...
            ,'TickLength',[.015 .015],'Layer','bottom')
        
        if e==1 && i==1
            plot([-1 -1],[mMin mMin+20],'k','LineWidth',3)
            plot([-1 0],[mMin*1.5 mMin*1.5],'k','LineWidth',3)
            
        end
        count=count+1;
        
    end
end


match_ylim(ax)
%set(gcf,'Position',[   680   248   109   730]);


subplot(1,7,7)
for e=1:3
    curve=squeeze(mean(BigExpect(N,TIO,:,e),2))*100;
    curveSE=squeeze(BigExpectation(N,:,e))*100;
    
    errorbar(Ori,curve,curveSE','o','Color',...
        cols(e,:),'MarkerFaceColor',cols(e,:),'CapSize',0,'LineWidth',2)
    
    hold on
    
    
    y=curve;
    min_y=min(y);
    max_y=max(y);
    spread=(max_y-min_y);
    [max1, ind1] = max(y);
    
    lb = [-100 0 0 10];
    ub = [1000 spread*1.15 180 100];
    x0=double([min_y  spread Ori(ind1) 30]);
    %ylim([-.010 .020])
    psy1 = lsqcurvefit(@fitCic180,x0,Ori,y',lb,ub);
    xticks(0:90:180)
    fitLine1 = fitCic180(psy1,X1);
    
    plot((X1),fitLine1,'-','Color',cols(e,:),'LineWidth',2)
    
    
end
xlabel('Orientation (\circ)')
ylabel('dF/F ( % )')

box off
l=legend('Random','','Expected','','Unexpected','');
set(l,'Position',[ 0.0278482138410421
    0.272613612194257
    0.110144927536231
    0.563186813186813]);
set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom','FontWeight','normal')
set(gcf,'Position',[         109         792        1035         182]);
MakePrettyFigure
cd(figPath)
print(['ExampleNeuronsNewCurve-' num2str(N) '.eps'],'-depsc')


%%
Ncount=1;
clear Shifted
for N=sigN
    datRandom = squeeze(mean(BigExpect(N,TIO,:,3),2));
    
    [a ind1]=max(datRandom);
    for expect =[ 1:5]
        
        t=squeeze(BigExpect(N,:,:,expect));
        Shifted(Ncount,expect,:,:)=circshift(t,4-ind1,2);
    end
    
    
    
    Ncount=Ncount+1;
end
%%

%%
close all
datM= squeeze(mean(Shifted));
datSE =squeeze(std(Shifted))/sqrt(size(Shifted,1));
dat1=squeeze(datM(1:3,:,:))*100;

mMin=min(min(min(dat1)));
mMax=max(max(max(dat1)));
dif =mMax-mMin;
%legend('R','expected','UnExpceted')
count=1;
clear ax
for i=1:6
    ax(i)= subplot(1,7,i);
    
    fh=fill([0.25 ,0.25 ,1,1],[mMax mMin mMin mMax],[.85 .85 .85],'EdgeColor','none');
    
    for e=[1:3]
        dat=squeeze(datM(e,:,:))*100;
        dS =squeeze(datSE(e,:,:))*100;
        
        hold on
        
        shadedErrorBar(NewTime,dat(:,i),dS(:,i),{'Color',cols(e,:)},.7)
        
        plot([NewTime(1) NewTime(end)],[0 0 ],'k','LineWidth',1)
        plot([0 0],[mMin mMax],'k','LineWidth',1)
        
        
        xlim([-1 2.0])
        
        
        axis off
        box off
        set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
            'YAxisLocation','left','LineWidth',1.5,...
            'XMinorTick','off','YMinorTick','off'...
            ,'TickLength',[.015 .015],'Layer','bottom')
        
        if e==1 && i==1
            plot([-1 -1],[mMin mMin+20],'k','LineWidth',3)
            plot([-1 0],[mMin*1.5 mMin*1.5],'k','LineWidth',3)
            
        end
        count=count+1;
        
    end
end


match_ylim(ax)
%set(gcf,'Position',[   680   248   109   730]);


subplot(1,7,7)
for e=1:3
    dat = squeeze(mean(Shifted(:,e,TIO,:),3))*100;
    datM= squeeze(mean(dat));
    datSE =squeeze(std(dat))/sqrt(size(dat,1));
    
    errorbar(Ori,datM,datSE','o','Color',...
        cols(e,:),'MarkerFaceColor',cols(e,:),'CapSize',0)
    
    hold on
    
    
    y=datM';
    min_y=min(y);
    max_y=max(y);
    spread=(max_y-min_y);
    [max1, ind1] = max(y);
    
    lb = [-100 0 0 10];
    ub = [1000 spread*1.15 180 100];
    x0=double([min_y  spread Ori(ind1) 30]);
    %ylim([-.010 .020])
    psy1 = lsqcurvefit(@fitCic180,x0,Ori,y',lb,ub);
    xticks(0:90:180)
    fitLine1 = fitCic180(psy1,X1);
    
    plot((X1),fitLine1,'-','Color',cols(e,:),'LineWidth',2)
    
    
end
xlabel('Orientation (\circ)')
ylabel('dF/F ( % )')

box off
l=legend('Random','','Expected','','Unexpected','');
set(l,'Position',[ 0.0278482138410421
    0.272613612194257
    0.110144927536231
    0.563186813186813]);
set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom','FontWeight','normal')
set(gcf,'Position',[         109         792        1035         182]);
set(gcf,'Renderer','painters')
MakePrettyFigure
cd(figPath)
print(['Shifted-' num2str(N) '.eps'],'-depsc')
%%
%%

Cols = linspecer(5);

dat = squeeze(mean(BigExpect(sigN,:,:,1:3),4));
neurons = 1:length(sigN);
dat=squeeze(mean(dat(:,TIO,:),2));
[R sorted]=max(dat,[],2);
close all
figure('Color','w')
clear ax
for expect =[ 3 1 2]
    
    ax(expect)=subplot(3,1,expect);
    dat = squeeze(BigExpect(sigN,:,:,expect))*100;
    datPreferred=zeros(size(dat,2),size(dat,2));
    for N=1:size(dat,1);
        
        datPreferred(N,:)=squeeze(dat(N,:,sorted(N)));
    end
    datPreferred=filtfast(datPreferred,2,[],'gaussian',1);
    if expect==3
        [a order]=sort(mean(datPreferred(:,TIO),2),'descend');
    end
    imagesc(NewTime,neurons,real((datPreferred(order,:))),[-10 50])
    colorbar
    hold on
    plot([ 0 0],[neurons(1) neurons(end)])
    xlabel('Time (s)')
    xlim([-1 3])
    box off
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
    colormap(linspecer(50))
end
match_clim(ax)
set(gcf,'Position',[   251   125   299   861])
cd(figPath)
print(gcf, 'PreferredOrientationREsponse.pdf')

%%
close all
cols = linspecer(3);
clear h
c=1;
for e=[1  3 ]
    dat=(squeeze(BigS(e,SigB,:)));
    datM=mean(dat);
    datSE=std(dat)/sqrt(size(dat,1));
    
    h(c)=shadedErrorBar(NewTime,datM,datSE,{'Color',cols(e,:)},.5)
    hold on
    c=c+1;
end
dat=(squeeze(BigS([1 3],SigB,:)));
[datobs, datrnd] = cluster_test_helper([squeeze(dat(1,:,:))-squeeze(dat(2,:,:))]', 2000);

[h_mem, p_mem, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);



pclustu = unique(p_mem);
npclust = nnz(pclustu < 0.05);
for ipclust = 1:npclust % extract time range of each significant cluster and show in figure
    currind  = p_mem == pclustu(ipclust);
    plot([min(NewTime(currind)),max(NewTime(currind))],[.1 .1 ],'Color','k','LineWidth',3)
    
end

yticks(.1:.05:2)
plot([0 0],[0.07 .22],'k','LineWidth',1)
c=legend([h.mainLine  ],{'Random','Unexpected'});
xlim([-.5 1.2])
ylabel('Circular mean')
xlabel('Time (s)')
box off
%%axis square
set(gcf,'Position',[   680   322   426   656])
set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
    'YAxisLocation','left','LineWidth',1.5,...
    'XMinorTick','off','YMinorTick','off'...
    ,'TickLength',[.015 .015],'Layer','bottom')
set(gcf,'Renderer','painters')
MakePrettyFigure
cd(figPath)
print(['Circular Mean Modulated Neurons' '.eps'],'-depsc')
%%
%%
Ncount=1;
clear fitStore diffStore
TIO = NewTime>.25 & NewTime<1;
Curve= collapse(BigExpect(:,TIO,:,:),2);
options = optimoptions('lsqcurvefit','Display','off');

for N=sigN
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
        
        lb = [-100 -1 0 10];
        ub = [10000 spread*1.15 180 100];
        x0=double([min_y  spread Ori(ind1) 30]);
        
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub,options);
        
        fitStore(Ncount,expect,:)=psy1;
        
    end
    ControlCounter(Ncount)=expect;
    Ncount=Ncount+1;
end
%%

jointMo=[];
for n=1:length(sigN)
    
    if sum(SigB==n)==1
        jointMo(n)=1;
    else
        jointMo(n)=0;
        
    end
end
%%
clear ax
for var=2
    close all
    hold on
    
    
    
    dat= fitStore(:,[1 3],var);
    
    if var ==2
        dat=dat*100;
        dat=log(dat);
    end
    for n=1:size(dat,1)
        if jointMo(n) ==1
            c='r';
        else
            c='k' ;
        end
        plot(dat(n,1),dat(n,2),'.','Color',c,...
            'MarkerFaceColor',c,'MarkerSize',12)
        hold on
    end
    hold on
    plot([-1.5 6],[-1.5 6],'k')
    
    axis square
    
    [H,P,CI,stats] =ttest(dat(:,1),dat(:,2));
    title(sprintf('N = %2.0f | t = %2.2f, p = %2.4f',...
        size(fitStore,1),real(stats.tstat),real(P)))
    box off
    real(stats.tstat),
    ylabel('Amplitude expected')
    
    ylabel('Amplitude unexpected')
    
    box off
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom')
end
if var ==2
    xlim([-1 6])
    ylim([-1 6])
    yticks(-2:2:6);
    xticks(-2:2:6);
end
% end
set(gcf,'Position',[   330   603   358   307])
cd(figPath)
%print(['SummaryPredictionStatsSpikes' '.pdf'],'-dpdf')
export_fig(['ScatterPrectvsRand-',num2str(var)])
%%
c=1;
close all
for var =[2 4]
    
    d=fitStore(:,[1 2 3],var);
    
    if var ==2
        d=d*100;
        % dat=log(dat);
    end
    datM=mean(d);
    [H,P,CI,stats] =ttest((squeeze(d(:,1))),(squeeze(d(:,end))))
    
    datSE =std(d)/sqrt(length(d));
    subplot(1,2,c)
    c=c+1;
    
    %  boxplot(real(d))
    if var ==2
        Bv=0;
    else
        Bv=0;
    end
    barwitherr(datSE,datM,'FaceColor',[.7 .7 .7],'basevalue',Bv)
    set(gca,'XTickLabel',{'Random','Expected','Unexpected'})
    %axis([.5 2.5 -1 0])
    xlim([.5 3.5])
    if var ==2
        %   ylim([-1.5 0])
    end
    if var ==2
        %   title('Gain')
        ylabel('Gain ( %dF/F )')
        %  yticks(0:.05:.15)
        % ylim([0 2.5])
    elseif var==4
        %     title('Width')
        ylabel('Width (\circ)')
        yticks(0:10:30)
        % ylim([0 .15])
    else
        
    end
    title(sprintf('N = %2.0f \n t = %2.2f, p = %2.4f',...
        size(dat,1),real(stats.tstat),real(P)))
    
    xtickangle(45)
    box off
    set(gca,'FontSize',14,'TickDir','in','XAxisLocation','bottom',...
        'YAxisLocation','left','LineWidth',1.5,'FontWeight','normal',...
        'XMinorTick','off','YMinorTick','off'...
        ,'TickLength',[.015 .015],'Layer','bottom','FontName','Helvetica')
end
%suptitle(sprintf('N=%2.0f',length(d)))
set(gcf,'Position',[   680   653   415   445])
cd(figPath)
MakePrettyFigure
%print(['SummaryPredictionStatsSpikes' '.pdf'],'-dpdf')
print(gcf, 'SummaryPredictionStats.pdf')
%%
var=2
close all
for compare =1:3
    subplot(1,3,compare)
    if compare ==1
    dat= fitStore(ControlCounter==5,[1 3],var);
    elseif compare ==2
            dat= fitStore(ControlCounter==5,[1 5],var);

    else
            dat= fitStore(ControlCounter==5,[2 4],var);

    end
    if var ==2
        dat=dat*100;
        dat=log(dat);
    end
    for n=1:size(dat,1)
    
            c='k' ;
      
        plot(dat(n,1),dat(n,2),'.','Color',c,...
            'MarkerFaceColor',c,'MarkerSize',12)
        hold on
    end
    hold on
    plot([-1.5 6],[-1.5 6],'k')
    
    axis square
    
    [H,P,CI,stats] =ttest(dat(:,1),dat(:,2));
    title(sprintf('N = %2.0f | t = %2.2f, p = %2.4f',...
        size(dat,1),real(stats.tstat),real(P)))
    box off
    real(stats.tstat),
    if compare ==1
   xlabel('Amplitude random')
    
    ylabel('Amplitude unexpected')
    
    elseif compare ==2
            xlabel('Amplitude random')
    
    ylabel('Amplitude unexpected control')
    else
        
           xlabel('Amplitude expected')
    
    ylabel('Amplitude expected control')
        
    end
    
    box off
MakePrettyFigure
xlim([-3 6])
ylim([-3 6])

end
set(gcf,'Position',[    0.0165    0.5470    1.1570    0.4190]*750)
cd(figPath)
MakePrettyFigure
%print(['SummaryPredictionStatsSpikes' '.pdf'],'-dpdf')
print(gcf, 'ControlStats.pdf')
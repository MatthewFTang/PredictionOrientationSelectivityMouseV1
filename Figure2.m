setUpPaths
%%
for len=1:23
   processData
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

        tempStore = Store(trials,:,:);
        for ori =1:6
            d=(tempStore(tempOris==ori,:,:));
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
save('SigPredictionCirc','StoreMod','StoreExpect','StoreStats','StoreCirc','StoreStats','StoreStatTime','StoreExpectation')
%%
cd(outputPath)
load('SigPredictionCirc')
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


SigB=find(BigMod<.025);
cd(outputPath

save('SigB','SigB')


fs =29.133;
NewTime = -1/fs*60:1/fs:1/fs*100;

TIO = NewTime >0.250 & NewTime <1.000;

sigN=find(any(BigStat([1 3],:)<.05));

TIO = NewTime>0.25 & NewTime <1;

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

%% Figure 2A
close all
[sortR sortI]=sort(BigStat(3,:)');
N=sortI(6)
X1=0:180;
Ori=0:30:150;
TIO = NewTime >.25 & NewTime <1;
cols = linspecer(3);
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
        dat_or = dat;
        win_size = 2;
        for time = 1:size(dat,1)-win_size
            dat(time,:)=mean(dat_or(time:time+win_size,:));
        end

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
saveas(['ExampleNeuronsNewCurve-' num2str(N) '.eps'],'-depsc')
%% find the preferred orientation in the random condition and re-align
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
%% Figure 2B
close all
dat_or = Shifted;
ShiftedM=Shifted;
win_size = 5;
for time = 1:size(ShiftedM,3)-win_size
    ShiftedM(:,:,time,:)=mean(dat_or(:,:,time:time+win_size,:),3);
end


datM= squeeze(mean(ShiftedM));
datSE =squeeze(std(ShiftedM))/sqrt(size(Shifted,1));
dat1=squeeze(datM(1:3,:,:))*100;

mMin=min(min(min(dat1)));
mMax=max(max(max(dat1)));
dif =mMax-mMin;
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

        temp_dat = squeeze(Shifted(:,e,:,i));
        [H,P,CI,STATS] = ttest(temp_dat);

        if e==1
            height =14;
        elseif e ==2
            height = 15;
        else
            height =16;
        end
        h_mem= P < 0.05/(6*3*161);
        if any(h_mem)
            plot([NewTime(h_mem)],ones(1,sum(h_mem))*height,'.','Color',cols(e,:))
        end


        count=count+1;

    end
end

match_ylim(ax)

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

%% Figure 2F
c=1;
close all
for var =[2 ]

    d=fitStore(:,[1 2 3],var);

    if var ==2
        d=d*100;
        % dat=log(dat);
    end
    datM=mean(d);
    subplot(1,2,c)
    plot(d(:,1),d(:,3),'o','Color','k','MarkerSize',4,'MarkerFaceColor','k  ')
    hold on
    min_d = min(min(d))
    max_d = max(max(d))
    plot([min_d,max_d],[min_d,max_d],'k','LineWidth',1)
    MakePrettyFigure
    if var ==2
        title('Gain ( %dF/F )')
        x = logspace(0,100)
        plot(x,x,'k','LineWidth',2)
        set(gca,'Xscale','log','YScale','log')
        yticks([1,10,100])
        xticks([1,10,100])

        xlim([0,100])
        ylim([0,100])
    elseif var==4
        %     title('Width')
        title('Width (\circ)')
        xticks(0:30:60)
        yticks(0:30:60)
        xlim([10,70])
        ylim([10,70])

    end
    xlabel('Random')
    ylabel('Unexpected')

    axis square
    c=c+1
end
a1=gca;
dat =d(:,1)-d(:,3);
edges = linspace(min(dat),-min(dat),30)
subplot(1,2,2)
n = histcounts(dat  ,edges );
bar(edges(1:end-1),n,'FaceColor' ...
    ,[.5,.5,.5])
hold on
plot([0,0],[0,120],'--k','Color','k')
MakePrettyFigure
axis square


cd(figPath)
saveas(gcf, 'Gain_and_width_scatter.pdf')
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
saveas(gcf, 'ControlStats.pdf')
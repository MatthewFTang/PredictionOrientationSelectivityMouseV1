clearvars
setUpPaths
dataPath=[dataPath '/Anaesthesia'];
cd(dataPath);
fileNames = dir('*-SaveEpochedNew.mat');
%%
for len=1:length(fileNames)
    processDataAnaesthesia
    %%
    fs=29.133;
    clear ExpectMeaned Sig fano ExpectSE
    clear pStore
    for expect =1:ExpCond
        NewTime = -1/fs*60:1/fs:1/fs*100;

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
        TIO = NewTime >0.250 & NewTime <1.000;

        tempStore = Store(trials,:,:);
        clear Expectataion
        for ori =1:6
            d=(tempStore(tempOris==ori,:,:));
            d1=mean(d(:,:,TIO),3);

            Expectataion(:,ori,expect)=std(d1)/sqrt(size(d,1));
            ExpectMeaned(:,:,ori,expect)=collapse(d,1);
            ExpectSE(:,:,ori,expect)=squeeze(std(d,1)/sqrt(size(d,1)));

        end
        %%
        tempStoreNoTime = mean(tempStore(:,:,TIO),3);
        for N=1:size(tempStore,2)
            [p h comp]= myanova(tempStoreNoTime(:,N),tempOris);
            pStore(expect,N,:)=p;
        end

    end
    cd(outputPath)
    save([fileNames(len).name(1:end-19) '-TuningSignifcanceTestingNewSortingAna'],'ExpectMeaned','pStore','Expectataion')
end
%%
cd(outputPath)
td=dir('*TuningSignifcanceTestingNewSortingAna.mat');

OverallExpectedC=[];
SigNeuronsOverallC=[];
OverallExpectedC=nan(1000,161,6,6);
OverallExpectedSE=nan(1000,161,6,6);

SigNeuronsOverallC=nan(5,1000);
FanoStore=[];
c=1;
for i=1:length(td)
    load(td(i).name)

    td(i).name
    N=size(ExpectMeaned,1);
    OverallExpectedC(c:c+N-1,:,:,1:size(ExpectMeaned,4))=ExpectMeaned;

    SigNeuronsOverallC(:,c:c+N-1)=pStore;

    SessStore(c:c+N-1,:)=i;
    c=c+N;
end
OverallExpectedC=OverallExpectedC(1:c,:,:,:);
SigNeuronsOverallC=SigNeuronsOverallC(:,1:c)
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
fs=1/29.99
NewTime = -60*fs:fs:100*fs;
SigNeuronsOverall=SigNeuronsOverallC<.05;

S2=(squeeze(SigNeuronsOverall([1 3],:,:)));
Sig=find((any(S2)))


Ncount=1;
clear fitStore diffStore
TIO = NewTime>0.25 & NewTime<1;
Curve= squeeze(mean(OverallExpectedC(:,TIO,:,:),2));
Ori=0:30:150

for N=Sig
    t=[1 2 3];

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

        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub);

        fitStore(Ncount,expect,:)=psy1;

    end
    ControlCounter(Ncount)=expect;
    Ncount=Ncount+1;
end
%%
%%
Ncount=1;
clear fitStoreAwake
TIO = NewTime>0.25 & NewTime<1;
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

        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub);

        fitStoreAwake(Ncount,expect,:)=psy1;

    end
    ControlCounter(Ncount)=expect;
    Ncount=Ncount+1;
end
Ncount


%% Figure 6E
close all
var=2;
Ana=log(fitStore(:,[1  3],var));
AwaSuprise =Ana(:,2)-Ana(:,1);


Awake=log(fitStoreAwake(:,[1  3],var));
AwakeSuprise =Awake(:,2)-Awake(:,1);
Combined_data = real([AwakeSuprise;AwaSuprise])
conditions = [zeros(1,length(AwakeSuprise)),ones(1,length(AwaSuprise))]
vs=violinplot(Combined_data,conditions,'MedianColor',[0,0,0],'ShowMean', ...
    true,'BoxColor',[0,0,0],'LineWidth',16,'ViolinColor',[.5,.5,.5])
datM= [mean(AwaSuprise) mean(AwakeSuprise)];
datSE=[std(AwaSuprise)/sqrt(length(AwaSuprise)) std(AwakeSuprise)/sqrt(length(AwakeSuprise))];
[H,P,CI,stats] =ttest2(real(AwaSuprise),real(AwakeSuprise))


set(gca,'XTickLabel',{'Awake','Anaesthetized'},'FontWeight','bold')
ylabel('Suprise effect','FontWeight','bold')


set(gcf,'Position',[   680   653   261   445])
title(sprintf('N = %2.0f | t = %2.2f, p = %2.4f',...
    size(fitStore,1),real(stats.tstat),real(P)))
cd(figPath)
MakePrettyFigure
set(gcf,'Renderer','painters')
MakePrettyFigure
cd(figPath)
saveas(gcf,['AwakeVsAna' '.eps'])
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
%% Figure 6A
close all
Oris=(0:30:150)-90;
clear h
cols = linspecer(6);
for o=1:6
    subplot(1,6,o)
    mMax=8;
    mMin=-3;
    fh=fill([0.25 ,0.25 ,1,1],[mMax mMin mMin mMax],[.85 .85 .85],'EdgeColor','none');
    hold on
    for e=[2 1 3]

        dat= squeeze(Shifted(:,e,:,o));
        datM= squeeze(nanmean(dat));
        datSE = squeeze(nanstd(dat))/sqrt(length(dat));
        h(e)=shadedErrorBar(NewTime,datM,datSE,{'Color',cols(e,:),'LineWidth',1},.5);
        xlim([-.7 2])
        ylim([-10 10])
        yticks(-10:10:20)

        hold on

    end
    if  o==1
        plot([-0 0],[0 5],'k','LineWidth',3)
        plot([-1 0],[-4 -4 ],'k','LineWidth',3)

    end
    title(Oris(o))

    plot([-20 20],[0 0],'k')
    plot([-.1 -.1],[-100 100],'k')

    box off
    ylim([-4 8])
    xticks([0 3])
    axis off
end

set(gcf,'Position',[  109.0000  792.0000  697.3333  182.0000]);
set(gcf,'Renderer','painters')
cd(figPath)
print(['AlignedNeuronsAna' '.eps'],'-depsc')
%% Figure 6B
X1=-90:90;
OriPreferred = -90:30:60;
close all
for e=1:3
    dat = squeeze(mean(Shifted(:,e,TIO,:),3));
    datM= squeeze(mean(dat));
    datSE =squeeze(std(dat))/sqrt(size(dat,1));

    errorbar(OriPreferred,datM,datSE','o','Color',...
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
    psy1 = lsqcurvefit(@fitCic180,x0,OriPreferred,y',lb,ub);
    xticks(OriPreferred)
    fitLine1 = fitCic180(psy1,X1);

    plot((X1),fitLine1,'-','Color',cols(e,:),'LineWidth',2)


end
xlabel('Orientation (\circ)')
ylabel('dF/F ( % )')

box off

set(gcf,'Position',[  109.0000  751.0000  220.6667  223.0000]);
set(gcf,'Renderer','painters')
MakePrettyFigure
cd(figPath)
print(['AnaTuningPreferred' '.eps'],'-depsc')

%% Figure 6f

Cols = linspecer(5);

close all
figure('Color','w')
clear ax
for expect =[3 1 2]

    ax(expect)=subplot(3,1,expect);
    dat = squeeze(Shifted(:,expect,:,4));
    neurons=1:size(dat,1);
    datPreferred=filtfast(dat,2,[],'gaussian',1);
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
set(gcf,'Renderer','painters')
print(['ResponseTimeCoursePreferredAna' '.eps'],'-depsc')
%%
%%
c=1;
close all
for var =[1:4]

    d=fitStore(:,[1 2 3],var);

    if var ==1

        d=real((d*100));
        t=  'baseline';
    elseif var ==2
        t= 'gain';
        d=real((d*100));

    elseif var ==3
        t= 'preferred orientation';
    else
        t= 'width';
    end
    datM=mean(d);
    [H,P,CI,stats] =ttest((squeeze(d(:,1))),(squeeze(d(:,end))))

    datSE =std(d)/sqrt(length(d));
    subplot(2,2,var)


    violinplot(d,{'Random','Expected','Unexpected'},'ViolinColor',cols)

    title(sprintf('%s \n t (%2.0f) = %2.2f, p = %2.4f',...
        t,stats.df,abs(stats.tstat),real(P)))
    if var ==1
        ylabel('dF/F ( % )')

    elseif  var ==2
        ylabel('dF/F ( % )')

        set(gca,'Yscale','log')
    elseif var ==3;
        ylabel('Orientation ( \circ )')
    else
        ylabel('Orientation ( \circ )')

    end
    xtickangle(45)
    MakePrettyFigure

end
set(gcf,'Renderer','painters')
set(gcf,'Position',[   614   415   666   723])
cd(figPath)
saveas(gcf, 'SummaryPredictionStatsAna.pdf')
%% Figure 6cd
close all
var =2
subplot(1,2,1)
    d=fitStore(:,[1 3],var);
      scatter(d(:,1),d(:,2),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')
      set(gca,'XScale','log','YScale','log')
      xlim([0.15,10000])
      ylim([0.15,10000])
      x =logspace(0,100);
      hold on
      plot(x,x,'k','LineWidth',2)
axis  square

xlabel('Random','FontWeight','bold')
ylabel('Unexpected','FontWeight','bold')
MakePrettyFigure
subplot(1,2,2)

histogram(d(:,1)-d(:,2),15,'FaceColor',[.5,.5,.5])
hold on 
plot([0,0],[0,40],'LineWidth',2)
xlim([-50,50])
axis square
MakePrettyFigure

cd(figPath)
saveas(gcf,'scatter_ana_gain.pdf')
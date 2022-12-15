clearvars
setUpPaths
%%
for len=1:23
   processData
    
    fs=29.133;
    
    NewTime = -1/fs*60:1/fs:1/fs*100;
    
    TIO = NewTime >0.250 & NewTime <1.000;
    baseTIO = NewTime>-.1 & NewTime<0;
    trials=StoreType==1;
    
    ind1 =trials & NewSupriseStore1==0;
    
    ind2= trials &NewSupriseStore1==1;
    storeNoTime = squeeze(mean(Store(:,:,TIO),3));
    baseTime = squeeze(mean(Store(:,:,baseTIO),3));
    response  = storeNoTime-baseTime;
    [H,P,CI,STATS]= ttest(response(trials&ind1,:));
    pStoreMismatch{len,1}=P;
    
    [H,P,CI,STATS]= ttest(response(trials&ind2,:));
    pStoreMismatch{len,2}=P;
    %%
    trials=StoreType==1 &NewSupriseStore1==1 ;
    tempOris = StoreX2(trials);
    expectedOri = ExpectedOrientationsActual(trials);
    uniOris = unique(expectedOri);
    tempStore = Store(trials,:,:);
    
    ExpectMeaned=zeros(size(tempStore,2),size(tempStore,3),length(uniOris),length(uniOris));
    
    for exOri = 1:length(uniOris)
        ind = expectedOri==uniOris(exOri);
        for ori =1:6
            d=(tempStore(tempOris==ori&ind,:,:));
            d= d-mean(d(:,:,1:60),3);
            d1=mean(d(:,:,TIO),3);
            ExpectMeaned(:,:,ori,exOri)=collapse(d,1);
            
        end
    end
    %%
    trials= RandTrials;
    tempOris = StoreX2(trials);
    tempStore = Store(trials,:,:);
    clear pStore
    tempStoreNoTime = mean(tempStore(:,:,TIO),3);
    for N=1:size(tempStore,2)
        [p, h, comp]= myanova(tempStoreNoTime(:,N),tempOris);
        pStore(N,:)=p;
    end
    
    storePstore{len}=pStore;
    storeExpectMeaned{len}=ExpectMeaned;
    
end
%%
cd(outputPath)
load('SigPredictionCirc')

BigStat=[];

N=1;
for len=1:length(StoreCirc)
    dat=StoreStats{len};
    newN=size(dat,2);

    BigStat(1:size(dat,1),N:N+newN-1)=dat;

    N=N+newN;
end
sigN=find(any(BigStat([1 3],:)<.05)); % which neurons show an orientation

%%
N=1

pBig = zeros(10000,2);
pStore=[];
BigExpectation=[]
for ii =1:length(pStoreMismatch)
    dat = pStoreMismatch{ii,1};
    nNew = length(dat);
    
    pBig(N:N+nNew-1,1)=dat;
    pBig(N:N+nNew-1,2)=pStoreMismatch{ii,2};
    BigExpectation(N:N+nNew-1,:,:,:)=storeExpectMeaned{ii};
    pStore(N:N+nNew-1)= storePstore{ii};
    N=N+nNew;
    
end

pBig=pBig(1:N-1,:);
%%
sig =pBig<.05;
sigError = sig(:,2)==1;
mean(sigError)
sum(sig(sigError,1))

%%
close all
% sigNeurons = pStore<.05;
sigNeurons = sigN;

expNoTime = (collapse(BigExpectation(sigNeurons,TIO,:,:),[2 4]));
[pre preOri ]=max(expNoTime,[],2);
uniPre=unique(preOri);
cols = linspecer(7,'green');
cols(1,:)=[];
order = (1:6)-3
cols = cols([1:2 6:-1:3],:);
Oris =0:30:150;
noTime =squeeze( mean(BigExpectation(sigNeurons,TIO,:,:),2));
x1= 0:180;
n=1;
for ii =1:length(uniPre);
    
    ind = preOri== uniPre(ii);
    dat= noTime(ind,:,:)*100;
    datM =squeeze(mean(dat,1));
    datSE= squeeze(std(dat)/sqrt(sum(ind)));
    ax(ii)= subplot(1,6,ii)
    hold on
    for eO =1:6
        
        y=datM(:,eO);
        bad= isnan(y);
        
        Ori = Oris (~bad);
        
        y=y(~bad);
        
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        
        lb = [-100 0 0 10];
        ub = [1000 spread*1.15 180 100];
        x0=double([min_y  spread Ori(ind1) 30]);
        
        errorbar(Oris,datM(:,eO),datSE(:,eO),'.','Color',cols(eO,:),'MarkerSize',24,'CapSize',0)
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y',lb,ub);
        xticks(0:90:180)
        fitLine1 = fitCic180(psy1,x1);
        
        plot(x1,fitLine1,'LineWidth',2,'Color',cols(eO,:))
        
        if eO ==3
            
            xlabel('Orientation')
            
            
        end
        if eO ==1;
            
            
            xlabel('dF/F ( % )')
        end
    end
%     title(sprintf('%2.0f | n = %2.0f',Oris(ii),))
    MakePrettyFigure
    if ii==1
        stringOri=num2str(uniOris');
        legend('30','','60','','90','','120','','150','','180')
    end
end
match_ylim(ax)

set(gcf,'Position',[    0.4383    0.9803    1.5440    0.3007]*1000)
cd(figPath)
print(['AllNeuronsdByExpected' '.eps'],'-depsc')
%%

%%
close all
expNoTime = (collapse(BigExpectation(sigNeurons,TIO,:,:),[2 4]));
[pre preOri ]=max(expNoTime,[],2);
uniPre=unique(preOri);
cols = linspecer(7,'green');
cols(1,:)=[];
order = (1:6)-3
cols = cols([1:2 6:-1:3],:);
Oris =0:30:150;
noTime =squeeze( mean(BigExpectation(sigNeurons,TIO,:,:),2));
x1= 0:180;
n=1;
for ii =1:length(uniPre);
    
    ind = preOri== uniPre(ii);
    dat= noTime(ind,:,:)*100;
    datM =squeeze(mean(dat,1));
    datSE= squeeze(std(dat)/sqrt(sum(ind)));
    ax(ii)= subplot(1,6,ii)
    hold on
    for eO =1:6
        
        y=datM(:,eO);
        bad= isnan(y);
        
        Ori = Oris (~bad);
        
        y=y(~bad);
        
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        
        lb = [-100 0 0 10];
        ub = [1000 spread*1.15 180 100];
        x0=double([min_y  spread Ori(ind1) 30]);
        
        errorbar(Oris,datM(:,eO),datSE(:,eO),'.','Color',cols(eO,:),'MarkerSize',24,'CapSize',0)
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y',lb,ub);
        xticks(0:90:180)
        fitLine1 = fitCic180(psy1,x1);
        
        plot(x1,fitLine1,'LineWidth',2,'Color',cols(eO,:))
        
        if eO ==3
            
            xlabel('Orientation')
            
            
        end
        if eO ==1
            
            
            xlabel('dF/F ( % )')
        end
    end
    title(sprintf('%2.0f \n n = %2.0f',Oris(ii),sum(ind)))
    MakePrettyFigure
    if ii==1
        stringOri=num2str(uniOris');
        legend('30','','60','','90','','120','','150','','180')
    end
end
match_ylim(ax)

set(gcf,'Position',[    0.4383    0.9803    1.5440    0.3007]*1000)
cd(figPath)
print(['AllNeuronsdByExpected' '.eps'],'-depsc')




%%
close all
% sigNeurons = pStore<.05;
expNoTime = (collapse(BigExpectation(sigNeurons,TIO,:,:),[2 4]));
[pre preOri ]=max(expNoTime,[],2);
uniPre=unique(preOri);
cols = linspecer(6);
Oris =0:30:150;
noTime =squeeze( mean(BigExpectation(sigNeurons,TIO,:,:),2));
x1= 0:180;
psyStore=[]
reAligned = [];
n=1;
realigned=zeros(size(noTime,1),6,6);

for ii =1:length(uniPre);
    
    ind = preOri== uniPre(ii);
    indFind= find(ind);
    dat= noTime(ind,:,:);
    dat=circshift(dat,4-ii,2);
    dat=circshift(dat,4-ii,3);
    subplot(1,6,ii)
    
    datM=squeeze(mean(dat));
    
    realigned(indFind,:,:)=dat;
    
end

realigned1(:,:,1)=realigned(:,:,3);
realigned1(:,:,3)=nanmean(realigned(:,:,[1 5]),3);
realigned1(:,:,2)=nanmean(realigned(:,:,[2 4]),3);
realigned1(:,:,4)=nanmean(realigned(:,:,[6]),3);
%%
opts = optimset('Display','off');

errorOri1=-60:30:90;
for N=1:size(realigned,1)
    disp(N)
    for e =1:6
        y=squeeze(realigned(N,:,e))*100;
        bad= isnan(y);
        
        Ori = errorOri1(~bad);
        
        y=y(~bad);
        
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        lb = [-100 0          -60  10];
        ub = [100 spread*1.15 90 100];
        x0=double([min_y  spread 0 30]);
        
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub,opts);
        psyStore(N,e,:)=psy1;
    end
end

%%
close all
clear exp
errorOri = Oris - 90;

cols1 = linspecer(8,'green');
cols1(1:2,:)=[];

cols(1,:)=cols1(2,:);
cols(2,:)=cols1(4,:);
cols(3,:)=cols1(6,:);
cols(4,:)=cols1(3,:);

cols(5,:)=cols1(2,:);
cols(6,:)=cols1(1,:);


Oris =0:30:150;
X=0:180
%cols =cols1

for e = 1:6
    
subplot(2,3,[1 4])
    
    ind = preOri== 4;
    dat= noTime(ind,:,:)*100;
    datM =squeeze(mean(dat,1));
    datSE= squeeze(std(dat)/sqrt(sum(ind)));
    hold on
    for eO =1:6
        
        y=datM(:,eO);
        bad= isnan(y);
        
        Ori = Oris (~bad);
        
        y=y(~bad);
        
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        
        lb = [-100 0 0 10];
        ub = [1000 spread*1.15 180 100];
        x0=double([min_y  spread Ori(ind1) 30]);
        
        errorbar(Oris,datM(:,eO),datSE(:,eO),'.','Color',cols(eO,:),'MarkerSize',24,'CapSize',0,'LineWidth',2)
        psy1 = lsqcurvefit(@fitCic180,x0,Ori,y',lb,ub);
        xticks(0:90:180)
        fitLine1 = fitCic180(psy1,X);
        
        plot(X,fitLine1,'LineWidth',2,'Color',cols(eO,:)) 
    end
    MakePrettyFigure
end
l1=legend('0','','30','','60','','90','','120','','150')
set(l1,'Position',[[0.2695 0.6647 0.0875 0.2992]])
x1=-90:60

%cols(5,:)=cols1(3,:);

for e =1:6
    subplot(2,3,[2 5])
    dat = realigned(:,:,e)*100;
    
    datM=squeeze(mean(dat));
    datSE=squeeze(std(dat)/sqrt(size(dat,1)));
    errorbar(errorOri,datM,datSE,'.','MarkerSize',24,'CapSize',0,'Color',cols(e,:),'LineWidth',2)
    hold on
    
    
    y=datM;
    bad= isnan(y);
    
    Ori = errorOri(~bad);
    
    y=y(~bad);
    
    min_y=min(y);
    max_y=max(y);
    spread=(max_y-min_y);
    [max1, ind1] = max(y);
    
    lb = [-100 0 0 10];
    ub = [1000 spread*1.15 180 100];
    x0=double([min_y  spread Ori(ind1) 30]);
    
    psy1 = lsqcurvefit(@fitCic180,x0,Ori,y,lb,ub)
    xticks(0:90:180)
    fitLine1 = fitCic180(psy1,x1);
    
    plot(x1,fitLine1,'LineWidth',2,'Color',cols(e,:))
    
    
    
    xticks(errorOri)
end
MakePrettyFigure
l2=legend('-60','','-30','','0','','30','','60','','90')
set(l2,'Position',[[0.4380 0.7026 0.0832 0.2992]])
xlabel('Preferred orientation')
ylabel('dF/F ( % )')

errorOri2=[-60,-30,0,30,60,90];
c=3;
for var=[2 1]
    subplot(2,3,c)
    c=c+3;
    dat=squeeze(psyStore(:,:,var));
    
    
    t =table(dat(:,1),dat(:,2),dat(:,3),dat(:,4),'VariableNames',{'1','2','3','4'});
    writetable(t,[num2str(var) '.csv'])
    dat=log(dat);
    datM=mean(dat,1);
    datM=exp(datM);
    datSE=std(dat,1)/sqrt(size(dat,1));
    datSE=exp(datSE);
    plot(errorOri2,datM,'-k','LineWidth',2)
    hold on
    for i =1:6
        errorbar(errorOri2(i),datM(i),datSE(i),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'CapSize',0,'MarkerSize',10,'LineWidth',2)
        
        
    end
    xlim([-90 91])
    xticks(errorOri2)
    axis square
    MakePrettyFigure
    if var ==1
        ylabel('Response anti-preferred (df/f)')
        xlabel('Difference with preferred orientation (deg)')
        
    else
        
        ylabel('Orientation gain (df/f)')
    end
end


set(gcf,'Position',[   123.6667  747.0000  910.0000  369.3333
])
MakePrettyFigure
cd(figPath)
print(['AllNeuronsAlignedByPreferredByExpected' '.eps'],'-depsc')




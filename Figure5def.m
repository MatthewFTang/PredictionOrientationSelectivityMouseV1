%%
%
clearvars
setUpPaths
for len=1:23
    
  processData
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
    for run =1:2
        if run ==1
            expect_ratio=0;
            adapt_ratio =.5;
        else
            expect_ratio=0.25;
            adapt_ratio =0.25;            
        end

        
        OngoingAdapation=ones(1,size(Channels_Form,1));
        
        for trial =1:length(ExpectedOrientationsActual)
            
            Test_Orientation=(OrientationsActual(trial));
            Expected_Orientation=(ExpectedOrientationsActual(trial));
            Channels_Form1=Channels_Form;
            Channels_Form1=Channels_Form1.*OngoingAdapation';
            Expected_Response=1-(Channels_Form1(:,Expected_Orientation))*expect_ratio;
            %
            Channels_Form1=Channels_Form1.*Expected_Response;
            
            
            Response_Form=Channels_Form1(:,Test_Orientation);
            OngoingAdapation= (1-(1-OngoingAdapation)/recoveryRatio);
            OngoingAdapation=OngoingAdapation.*abs(Response_Form.*adapt_ratio-1')';
            
            Response_FormStore(run,trial,:)=Response_Form;
            
        end
    end
    


close all
%%
Response_FormStoreRotated = zeros(size(Response_FormStore));
uniOri = unique(OrientationsActual);
for o =1:6
    ind =OrientationsActual==uniOri(o);
    dat = Response_FormStore(:,ind,:);
    Response_FormStoreRotated(:,ind,:)=circshift(dat,3-o,3);
end

%%
cols = linspecer(3);
for run =1:2
    for expect =1:3
        fs=29.133;
        
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
        
        dat = squeeze(mean(Response_FormStoreRotated(run,trials,:),2));
        hold on
        plot(dat,'o','Color',cols(expect,:),'Color',cols(expect,:))
        
        
        modelResultsStore(len,run,expect,:)=dat;
    end
end
end

%%
psyStore=[];
chanOris1=chanOris-90;

for l=1:23
    
    for e=1:3
         if e ==1
        run =1;
    else
        run =2;
    end
        dat = squeeze(modelResultsStore(l,run,e,:));
        
        y=dat;
        min_y=min(y);
        max_y=max(y);
        spread=(max_y-min_y);
        [max1, ind1] = max(y);
        
        lb = [-100 0 0 10];
        ub = [1000 spread*1.15 180 100];
        x0=double([min_y  spread chanOris1(ind1) 30]);
        psy1 = lsqcurvefit(@fitCic180,x0,chanOris1,y',lb,ub);
        
        psyStore(l,e,:)=psy1;
    end
    
end
%%



close all
X1=-90:90;
subplot(1,3,1)
for e=1:3
    if e ==1
        run =1;
    else
        run =2;
    end
    dat = squeeze(modelResultsStore(:,run,e,:));
    datM=mean(dat);
    datSE=std(dat)/sqrt(23);
    hold on
    errorbar(chanOris1,datM,datSE,'o','Color',cols(e,:),'MarkerFaceColor',cols(e,:))
    
    
    y=datM';
    min_y=min(y);
    max_y=max(y);
    spread=(max_y-min_y);
    [max1, ind1] = max(y);
    
    lb = [-100 0 0 10];
    ub = [1000 spread*1.15 180 100];
    x0=double([min_y  spread chanOris1(ind1) 30]);
    psy1 = lsqcurvefit(@fitCic180,x0,chanOris1,y',lb,ub);
    
    fitLine1 = fitCic180(psy1,X1);
    

    plot((X1),fitLine1,'-','Color',cols(e,:),'LineWidth',2)
    
    
end
conditions={'Random','Expected','Unexpected'};
xticks([-90:90:90])
xlabel('Preferred orientation')
ylabel('Response (a.u.)','FontWeight','bold')
MakePrettyFigure
c=2;
for var=[2 4]
    subplot(1,3,c)

    dat = psyStore(:,:,var);

      
    datM=mean(dat);
    datSE=std(dat)/sqrt(23);
    for i=1:3
        fprintf('var :%2.0f i = %2.0f | m = %2.2f sD = %2.2f\n ',var,i,datM(i),std(dat(:,i)))
%         bar(i,datM(i),'FaceColor',cols(i,:),'EdgeColor',cols(i,:))
%         hold on

        
%         errorbar(i,datM(i),datSE(i),'k','LineWidth',2,'CapSize',1)
    end
    violinplot(dat,conditions,'ViolinColor',cols,'ShowBox',false)

    xticks(1:3)
    if var ==2
        ylabel('Gain (a.u.)','FontWeight','bold')
        ylim([0.5,.7])
    else
        ylabel('Width ( deg )','FontWeight','bold')
                ylim([20,35])

    end
    xticklabels(conditions);
    MakePrettyFigure
    c=c+1;
end
set(gcf,'Position',[  450.3333  836.3333  777.3334  382.6667])
cd(figPath)
set(gcf,'Render','Painters')
saveas(gcf,['ModelResults.pdf'])

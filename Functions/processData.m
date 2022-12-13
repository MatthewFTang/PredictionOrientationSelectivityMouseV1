
%%
cd(dataPath);
load(fileNames(len).name)
fprintf('loading session : %2.0f\n',len)

%%
clear epochedFF
TrialCount=1;
epochedFF=nan(length(StoreType),sum(Iscell),161);
for exp =1:length(StoreTimes)
    Ff=StoreFF{exp};
    TrialTimes = StoreTimes{exp};
    CaTimeStamp=StoreCaTimes{exp};
    for trial=1:length(TrialTimes)
        ind=find(CaTimeStamp'>TrialTimes(trial),1);

        temp =Ff(ind-60:ind+100,logical(Iscell))';

        epochedFF(TrialCount,:,:)=temp;

        TrialCount=TrialCount+1;
    end

end
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
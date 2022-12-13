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
    if Nexps ==4
        ExpCond=3;

    else
        ExpCond=5;
    end
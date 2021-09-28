%% finds the suprising data in the same way as the experimental analysis
%close all
oriSequence = (OrientationsStore);

diffOri=[1 diff(oriSequence)];
modDiffOri=mod(diffOri,6);

% plot(oriSequence,'-o')
% hold on
% plot(modDiffOri,'-o')

clear ExpectedA
ExpectedA(1:2)=1;
ExpectedOrientation =[];

%ExpectedOrientation(ExpectedOrientation==8)=2;

%ExpectedOrientation(ExpectedOrientation==-1)=5;

NewSupriseStore =[0 diff(modDiffOri)~=0];
RunCounter=[];
NewSupriseStore1=[];
run=0;
for i=1:length(NewSupriseStore)
    if NewSupriseStore(i) ==1 && NewSupriseStore(i-1)~=1
        
        NewSupriseStore1(i)=1;
        run=0;
    else
        NewSupriseStore1(i)=0;
        run=run+1;
        
    end
    RunCounter(i)=run;
end
%plot(NewSupriseStore1,'g')
%plot(NewSupriseStore,'r')

for trial=2:length(oriSequence)
    if NewSupriseStore1(trial)==1
        if modDiffOri(trial-1)==modDiffOri(trial)
            ExpectedA(trial)=1;
            ExpectedOrientation(trial)=oriSequence(trial);
        else
            ExpectedA(trial)=0;
            if modDiffOri(trial-1)==5
                ExpectedOrientation(trial)=oriSequence(trial-1)-1;
            else
                ExpectedOrientation(trial)=oriSequence(trial-1)+1;
                
            end
        end
    else
        ExpectedOrientation(trial)=oriSequence(trial);
        
    end
    if StoreType(trial)==2
        ExpectedOrientation(trial)=randperm(6,1);
    end
end

ExpectedOrientation(ExpectedOrientation==7)=1;
ExpectedOrientation(ExpectedOrientation==0)=6;
%plot(RunCounter,'-ko')
%plot(StoreType)
%plot(ExpectedOrientation,'-ro')
%xlim([7200 7300])
TestOri=0:30:150;
ExpectedOrientationsActual = TestOri(ExpectedOrientation);
ExpectedOrientationsActual(ExpectedOrientationsActual==0)=180;
%set(gcf,'Position',[   236   565   950   240])
%%
clockwiseSeq = [1:6 1:6];
anticlockwiseSeq=[6:-1:1 6:-1:1];
% close all
% plot(oriSequence,'-go')
% hold on
%   plot(diffOri,'-bo')

% plot(modDiffOri,'-ro')
%xlim([7200 7300])
clear cat
distance=nan(size(oriSequence));
for trial =5:length(oriSequence)
    if (oriSequence(trial-1)-oriSequence(trial-2))==1 ...
            ||(oriSequence(trial-1)-oriSequence(trial-2))==-5
        a=clockwiseSeq(oriSequence(trial-1):end);
        distance(trial)=  find(oriSequence(trial)==a,1)-1;
        
    elseif (oriSequence(trial-1)-oriSequence(trial-2))==-1 ...
            ||(oriSequence(trial-1)-oriSequence(trial-2))==5
        a=anticlockwiseSeq(find(oriSequence(trial-1)==anticlockwiseSeq,1):end);
        distance(trial)=  find(oriSequence(trial)==a,1)-1;
    end
    
    
end
%   plot(distance,'o-')

TrialCount=1;
correctedDistance=distance;
for trial=1:length(distance)
    if correctedDistance(trial)~=1 && ~isnan(correctedDistance(trial))
        correctedDistance(trial+1:trial+4)=NaN;
    end
    
end
%plot(correctedDistance,'o-')

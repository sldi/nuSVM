% clear all
% for f=1:25
%     clearvars -except f Pd Pf Pdv Pfv par
% expl1= load('HMTD250 ng31data_DDT.mat_all_featurevector.mat', 'all_feature');
% expl2= load('HMTD50 ng33data_DDT.mat_all_featurevector.mat', 'all_feature');
% expl3= load('HMTD25 ng10data_DDT.mat_all_featurevector.mat', 'all_feature');
% expl4= load('HMTD10 ng31data_DDT.mat_all_featurevector.mat', 'all_feature');
% expls1= load('HMTD250 ng40data_SST.mat_all_featurevector.mat', 'all_feature');
% expls2= load('HMTD100 ng39data_SST.mat_all_featurevector.mat', 'all_feature');
% expls3= load('HMTD25 ng40data_SST.mat_all_featurevector.mat', 'all_feature');
% blank_train= load('INDCheckPoint0 ng1807data_IND_Operational Blanks.mat_all_featurevector.mat', 'all_feature');
% blank = load('INDCheckPoint0 ng1062data_IND2_Operational Blanks.mat_all_featurevector.mat', 'all_feature');
% expl1 = expl1.all_feature;
% expl2 = expl2.all_feature;
% expl3 = expl3.all_feature;
% expl4 = expl4.all_feature;
% expls1 = expls1.all_feature;
% expls2 = expls2.all_feature;
% expls3 = expls3.all_feature;
% blank_train = blank_train.all_feature;
% blank = blank.all_feature;
% HMTD_train = horzcat(expl1,expl2,expl3,expl4);
% HMTD = horzcat(expls1,expls2,expls3);
% 
% blank=blank';
% blank_train=blank_train';
% HMTD=HMTD';
% HMTD_train=HMTD_train';
% Xtrain = vertcat(HMTD_train,blank_train);
% Ytrain = vertcat(ones(size(HMTD_train,1),1),-1*ones(size(blank_train,1),1));
% 
% val_percentage=0.6;
% nr=size(HMTD,1);
% randindex=rand(nr,1);
% HMTD_test=HMTD(find(randindex>val_percentage),:);
% HMTD_val=HMTD(find(randindex<=val_percentage),:);
% nr=size(blank,1);
% randindex=rand(nr,1);
% blank_test=blank(find(randindex>val_percentage),:);
% blank_val=blank(find(randindex<=val_percentage),:);
clear all
for f=1:25
    clearvars -except f Pd Pf Pdv Pfv par
expl1= load('TATP250 ng11data_DDT.mat_all_featurevector.mat', 'all_feature');
expl2= load('TATP100 ng30data_DDT.mat_all_featurevector.mat', 'all_feature');
expl3= load('TATP25 ng20data_DDT.mat_all_featurevector.mat', 'all_feature');
expl4= load('TATP50 ng50data_DDT.mat_all_featurevector.mat', 'all_feature');
expls1= load('TATP100 ng20data_SST.mat_all_featurevector.mat', 'all_feature');
expls2= load('TATP250 ng20data_SST.mat_all_featurevector.mat', 'all_feature');
expls3= load('TATP50 ng20data_SST.mat_all_featurevector.mat', 'all_feature');
% expls4= load('TNT2.5 ng40data_SST.mat_all_featurevector.mat', 'all_feature');
% expls5= load('TNT25 ng40data_SST.mat_all_featurevector.mat', 'all_feature');
blank_train= load('INDCheckPoint0 ng1807data_IND_Operational Blanks.mat_all_featurevector.mat', 'all_feature');
blank = load('INDCheckPoint0 ng1062data_IND2_Operational Blanks.mat_all_featurevector.mat', 'all_feature');
expl1 = expl1.all_feature;
expl2 = expl2.all_feature;
expl3 = expl3.all_feature;
expl4 = expl4.all_feature;
expls1 = expls1.all_feature;
expls2 = expls2.all_feature;
expls3 = expls3.all_feature;
% expls4 = expls4.all_feature;
% expls5 = expls5.all_feature;
blank_train = blank_train.all_feature;
blank = blank.all_feature;
HMTD_train = horzcat(expl1,expl2,expl3,expl4);
HMTD = horzcat(expls1,expls2);
tester=expls1';
blank=vertcat(blank',blank_train');
HMTD=HMTD';
HMTD_train=HMTD_train';

val_percentage=0.6;
nr=size(HMTD,1);
randindex=rand(nr,1);
HMTD_tv=HMTD(find(randindex>val_percentage),:);
HMTD_con=HMTD(find(randindex<=val_percentage),:);
testindex=rand(size(tester,1),1);
testintermediate= tester(find(testindex>val_percentage),:);
testrain = tester(find(testindex<=val_percentage),:);
nr=size(blank,1);
randindex=rand(nr,1);
blank_tv=blank(find(randindex>val_percentage),:);
blank_con=blank(find(randindex<=val_percentage),:);


test_percentage=0.5;
nr=size(HMTD_tv,1);
randindex=rand(nr,1);
testindex=rand(size(testintermediate,1),1);
HMTD_test=testintermediate(find(testindex>test_percentage),:);
HMTD_val=HMTD_tv(find(randindex<=test_percentage),:);
testval=testintermediate(find(testindex<=test_percentage),:);

nr=size(blank_tv,1);
randindex=rand(nr,1);
blank_test=blank_tv(find(randindex>test_percentage),:);
blank_val=blank_tv(find(randindex<=test_percentage),:);


Xtrain = vertcat(HMTD_train,HMTD_con,testrain,blank_con);
Ytrain = vertcat(ones(size(HMTD_train,1),1),ones(size(HMTD_con,1),1),ones(size(testrain,1),1),-1*ones(size(blank_con,1),1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xval = vertcat(HMTD_val,blank_val);
Xtest = vertcat(HMTD_test,blank_test);
Yval = vertcat(ones(size(HMTD_val,1),1),-1*ones(size(blank_val,1),1));
Ytest = vertcat(ones(size(HMTD_test,1),1),-1*ones(size(blank_test,1),1));
d=size(Xtrain,2);
n=size(Xtrain,1);
np=length(find(Ytrain==1));
nn=n-np;

for k = 1:5;
for l = 1:5;
    nun=k/5;
    nup=l/5;
nu(k,l) = (2*nup*nun*nn*np)/(n*(nun*nn+nup*np));
gamma(k,l) = (nu(k,l)*n)/(2*nup*np);


for i=1:n
    if(Ytrain(i)==1)
        weighted_gamma(i)=gamma(k,l)/n;
    elseif(Ytrain(i)==-1)
        weighted_gamma(i)=(1-gamma(k,l))/n;
    end
end

%%%%%%%%%%%%%%%%%%%CVX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
%2nuSVM
    variables wtrain(d) e(n) btrain rho 
    dual variable alphatrain
    minimize( 0.5*wtrain'*wtrain - nu(k,l)*rho+sum(weighted_gamma'.*e))
    subject to
        Ytrain.*(Xtrain*wtrain+btrain)- rho +e >=0   :alphatrain;
        e>=0; %slack
        rho>=0;
% %CSVM
%     variables wtrain(d) e(n) btrain 
%     dual variable alphatrain
%     minimize( 0.5*wtrain'*wtrain + C*sum(e)) %norm(w) almost works except it takes an extra sqrt
%     subject to
%         Yt.*(Xt*wtrain+btrain)- 1 +e >=0   :alphatrain;
%         e>=0; %slack
cvx_end


%Primal Classification
 out = sign(wtrain'*Xval'+btrain);
%Dual Classification
% for i=1:n
%      out(i)=sum(alphatrain.*Yt.*(Xtest(i,:)*Xt')');
% end
% out = sign(out+btrain);
Pf_tempv=zeros(size(out));
Pd_tempv=zeros(size(out));
Pm_tempv=zeros(size(out));
%Assume 1 is explosive, -1 is blank
for i=1:length(out)
    if (out(i)==1 && Yval(i)==-1)
        Pf_tempv(i)=1;
    elseif (out(i)==1 && Yval(i)==1)
        Pd_tempv(i)=1;
    elseif (out(i)==-1 && Yval(i)==1)
        Pm_tempv(i)=1;
    end
end
Pf_tempv=sum(Pf_tempv)/length(find(Yval==-1));
Pm_tempv=sum(Pm_tempv)/length(find(Yval==1));
Pd_tempv=sum(Pd_tempv)/length(find(Yval==1));
Pdv(f)=Pd_tempv;
Pfv(f)=Pf_tempv;
ECV(k,l,1)=Pd_tempv;
ECV(k,l,2)=Pf_tempv;

end

end

max=0;
minPf = 0.003;
indicator=0;
for k = 1:5;
for l = 1:5;
if ECV(k,l,2)<=minPf && ECV(k,l,1)>=max
    indicator = 1;
    max=ECV(k,l,1);
    minPf=ECV(k,l,2);
    par{f} = [nu(k,l),gamma(k,l)];
end
end
end
if indicator==0
   continue
end
%%TEST TEST TEST TEST
for i=1:size(Xtrain,1)
    if(Ytrain(i)==1)
        weighted_gamma(i)=par{f}(1,2)/n;
    elseif(Ytrain(i)==-1)
        weighted_gamma(i)=(1-par{f}(1,2))/n;
    end
end
%%%%%%%%%%%%%%%%%%%CVX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
%2nuSVM
    variables wtrain(d) e(n) btrain rho 
    dual variable alphatrain
    minimize( 0.5*wtrain'*wtrain - par{f}(1,1)*rho+sum(weighted_gamma'.*e))
    subject to
        Ytrain.*(Xtrain*wtrain+btrain)- rho +e >=0   :alphatrain;
        e>=0; %slack
        rho>=0;
% %CSVM
%     variables wtrain(d) e(n) btrain 
%     dual variable alphatrain
%     minimize( 0.5*wtrain'*wtrain + C*sum(e)) %norm(w) almost works except it takes an extra sqrt
%     subject to
%         Yt.*(Xt*wtrain+btrain)- 1 +e >=0   :alphatrain;
%         e>=0; %slack
cvx_end


%Primal Classification
 out = sign(wtrain'*Xtest'+btrain);
%Dual Classification
% for i=1:n
%      out(i)=sum(alphatrain.*Yt.*(Xtest(i,:)*Xt')');
% end
% out = sign(out+btrain);
Pf_temp=zeros(size(out));
Pd_temp=zeros(size(out));
Pm_temp=zeros(size(out));
%Assume 1 is explosive, -1 is blank
for i=1:length(out)
    if (out(i)==1 && Ytest(i)==-1)
        Pf_temp(i)=1;
    elseif (out(i)==1 && Ytest(i)==1)
        Pd_temp(i)=1;
    elseif (out(i)==-1 && Ytest(i)==1)
        Pm_temp(i)=1;
    end
end
Pf_temp=sum(Pf_temp)/length(find(Ytest==-1));
Pm_temp=sum(Pm_temp)/length(find(Ytest==1));
Pd_temp=sum(Pd_temp)/length(find(Ytest==1));
Pf(f)=Pf_temp;
Pd(f)=Pd_temp;
end
figure
plot(Pd)
hold on
plot(Pf)
hold off

figure
scatter(Pf,Pd)
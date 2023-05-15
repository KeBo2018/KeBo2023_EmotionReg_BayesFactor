clear
load('D:\CANlab_Working\BF_Mat\Combined_001_3mm_Multiply_BF_4_6_WithoutRest3_Beta_WithoutSTS_NewBase_CC15')
load('D:\CANlab_Working\Behavioral variable Chart\RegNeg_Contrast_AHAB.mat')
RegNeg_AHAB=Reg_Neg;
RegNeg_AHAB=Reg_Neg;

% load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\AHAB_WithoutRest3_Whole_Beta.mat')
% % Whole_Neg_AHAB=Whole_Neg;
% Whole_Neg_AHAB=image_math(Whole_Neg,Whole_Neu,'minus');
% load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\PIP_WithoutRest3_Whole_Beta.mat')
% % Whole_Neg_PIP=Whole_Neg;

% Whole_Neg_PIP=image_math(Whole_Neg,Whole_Neu,'minus');


load('D:\CANlab_Working\Behavioral variable Chart\RegNeg_Contrast_PIP.mat')
RegNeg_PIP=Reg_Neg;


load('D:\CANlab_Working\Behavioral variable Chart\AHAB_rating.mat')
Reg_rating_AHAB=nanmean(Reg_Rating,2);
Neg_rating_AHAB=nanmean(Neg_Rating,2);
Neu_rating_AHAB=nanmean(Neu_Rating,2);
Success_AHAB=Neg_rating_AHAB-Reg_rating_AHAB;
EmoAct_AHAB=Neg_rating_AHAB-Neu_rating_AHAB;


% Success_AHAB=Neg_rating-Neu_rating;
% Success_AHAB=Neg_rating;
load('D:\CANlab_Working\Behavioral variable Chart\PIP_rating.mat')
Reg_rating_PIP=nanmean(Reg_Rating,2);
Neg_rating_PIP=nanmean(Neg_Rating,2);
Neu_rating_PIP=nanmean(Neu_Rating,2);
Success_PIP=Neg_rating_PIP-Reg_rating_PIP;
EmoAct_PIP=Neg_rating_PIP-Neu_rating_PIP;

data = [Neu_rating_PIP,Neg_rating_PIP,Reg_rating_PIP; Neu_rating_AHAB,Neg_rating_AHAB,Reg_rating_AHAB];
groups = repmat({'Reg '; 'Neg '; 'Neu'}, 100, 1);

colorcoding=[110 203 99;0 176 240; 255 89 76];
colorcoding=colorcoding/255;
% Create violin plot
figure
violinplot(data,'xlabel',{'Neu '; 'Neg '; 'Reg'},'facecolor',colorcoding,'edgecolor','k','bw',0.1,'pointsize',4, 'mc','k','facealpha',0.8,'linewidth',2, 'showMean', false, 'showMedian', false)
%               'medc','r--')
box off
set(gca,'fontsize',25,'fontweight','bold','LineWidth',3)
hold on
plot([1,2],data(:,1:2)','color',colorcoding(1,:))

hold on
plot([2,3],data(:,2:3)','color',colorcoding(3,:))

EmoAct_PIP=Neg_rating_PIP-Neu_rating_PIP;

% clear y x
% len = 200; sub = 20;
% x = zeros(len,sub);
% x(11:20,:) = 2;                   % create signal
% x(111:120,:) = 2;
% 
% x = mat2cell(x, size(x, 1), ones(1, size(x, 2)));
% y = mat2cell(y, size(y, 1), ones(1, size(y, 2)));
% 
% c = normrnd(0.5,0.1,sub,1);       % slope between-subjects variations
% d = normrnd(3,0.2,sub,1);         % intercept between-subjects variations
% 
% for i=1:sub, y(:,i) = d(i) + c(i).*x(:,i) + 2*i + normrnd(0,0.5,len,1); end
% 
% out = igls_multicond(y, x)  % for igls
% disp('Input random-effect variances: '); disp(std([d c]))
% disp('Est.  random-effect variances: '); disp(sqrt(out.betastar)');
% 
% 
% igls_multicond
Success_PIP=Neg_rating_PIP-Reg_rating_PIP;
% Success_PIP=Neg_rating-Neu_rating;
% Success_PIP=Neg_rating; 
Reg_Neg=image_math(RegNeg_AHAB,RegNeg_PIP,'concatenate')
% Whole_Neg=image_math(Whole_Neg_AHAB,Whole_Neg_PIP,'concatenate')
% Reg_Neg=image_math(RegNeg_AHAB,RegNeg_PIP,'concatenate')
Success=[Success_AHAB;Success_PIP];


Success_percentage= length(find(Success>0))/length(Success);

[h p ci tstat]=ttest(Success_PIP)
d=mean(Success_PIP)/std(Success_PIP);
[h p ci tstat]=ttest(Success_AHAB)
d=mean(Success_AHAB)/std(Success_AHAB);

[h p ci tstat]=ttest(EmoAct_PIP)
d=mean(EmoAct_PIP)/std(EmoAct_PIP);
[h p ci tstat]=ttest(EmoAct_AHAB)
d=mean(EmoAct_AHAB)/std(EmoAct_AHAB);

Success_Combined=[Success_AHAB;Success_PIP];
EmoAct_Combined=[ EmoAct_AHAB;EmoAct_PIP];

[h p ci tstat]=ttest(Success_Combined)
d=mean(Success_Combined)/std(Success_Combined);
[h p ci tstat]=ttest(EmoAct_Combined)
d=mean(EmoAct_Combined)/std(EmoAct_Combined);

% AcrossStudy=[ones(1,182), -1*ones(1,176)];
Brainmat=[mean(Reg_Neg.dat(indexReappraisalOnly,:),1); mean(Reg_Neg.dat(indexReappraisal_D,:),1)]
Brainmat1=[mean(Reg_Neg.dat(indexReappraisalOnly,:),1); mean(Reg_Neg.dat(indexReappraisal_D,:),1) ];
Brainmat2=[mean(Reg_Neg.dat(indexReappraisalOnly,:),1); ];
Brainmat3=[mean(Reg_Neg.dat(indexReappraisal_D,:),1); ];
% Brainmat=[ Brainmat ; AcrossStudy];

mdl=fitlm(Brainmat',Success);
corr(mdl3.Fitted,Success )

mdl1 = fitglm(Brainmat1(:,:)',Success)

mdl2 = fitlm(Brainmat2(:,:)',Success)

mdl3 = fitlm(Brainmat3(:,:)',Success)

RegNeg_AHAB.X=Success_AHAB;
RegNeg_PIP.X=Success_PIP;
Reg_Neg=image_math(RegNeg_AHAB,RegNeg_PIP,'concatenate')

% AHAB_Model=xval_SVR(RegNeg_AHAB.dat(indexReappraisal_D,:)',RegNeg_AHAB.X,(1:182)','norepeats', 'nobootstrap')
% PIP_Model=xval_SVR(RegNeg_PIP.dat(indexReappraisal_D,:)',RegNeg_PIP.X,(1:176)','norepeats', 'nobootstrap')

fitlinear(Brainmat1(:,:)',Success)
aic(mdl1)
[b,bint,r,rint,stats]=regress(Brainmat2(:,:)',Success)
corr(mean(Reg_Neg.dat(indexReappraisalOnly,:),1)', mean(Reg_Neg.dat(indexReappraisal_D,:),1)')
%mediation
load('D:\CANlab_Working\BF_Mat\AHAB_BF_4_6_WithoutRest3_Beta')
[R P]=corr(mean(RegNeg_AHAB.dat(indexReappraisal_D,:),1)', Success_AHAB)

M=Brainmat3;
X=Brainmat2;
Y=Success;
mdl1 = fitlm(X',Y') 
mdl2 = fitlm(X',M')
mdl3 = fitlm([X;M]',Y')
mdl4 = fitlm(M',Y')

[paths, stats] = mediation(X', Y, M')
%%%%%%%%%%Univariate%%%%%
[R(1) P(1)]=corr(mean(Reg_Neg.dat(indexOverlap,:),1)',Success);
[R(2) P(2)]=corr(mean(Reg_Neg.dat(indexReappraisalOnly,:),1)',Success);
[R(3) P(3)]=corr(mean(Reg_Neg.dat(indexReappraisal_D,:),1)',Success);
[R(4) P(4)]=corr(mean(Reg_Neg.dat(indexLookOnly,:),1)',Success);

for i=1:1000
selectingindex=randi(358,358,1);
[R P]=corr(mean(Reg_Neg.dat(indexOverlap,selectingindex),1)',Success(selectingindex));
RBoot(i)=R;
end
figure
histogram(RBoot,'BinWidth',0.01)

for i=1:1000
selectingindex=randi(358,358,1);
[R P]=corr(mean(Reg_Neg.dat(indexReappraisalOnly,selectingindex),1)',Success(selectingindex));
RBoot(i)=R;
end
figure
histogram(RBoot,'BinWidth',0.01)


for i=1:1000
selectingindex=randi(358,358,1);
[R P]=corr(mean(Reg_Neg.dat(indexReappraisal_D,selectingindex),1)',Success(selectingindex));
RBoot(i)=R;
end
figure
histogram(RBoot,'BinWidth',0.01)
% [R(1) P(1)]=corr(mean(Whole_Neg.dat(indexOverlap,:),1)',Success);
% [R(2) P(2)]=corr(mean(Whole_Neg.dat(indexReappraisalOnly,:),1)',Success);
% [R(3) P(3)]=corr(mean(Whole_Neg.dat(indexReappraisal_D,:),1)',Success);
% [R(4) P(4)]=corr(mean(Whole_Neg.dat(indexLookOnly,:),1)',Success);

%%%%%%%
% [R(1) P(1)]=corr(mean(Reg_Neg.dat(indexReappraisal_D,:),1)',Success,'type','spearman');

figure
dotcolor=[0.3010 0.7450 0.9330];
dotsize=500;
linewidth=3;
fontsize=12;
subplot(1,4,2)
scatter(mean(Reg_Neg.dat(indexOverlap,:),1)',Success,dotsize,dotcolor,'.');
h=lsline

set(h(1),'color','#FF594C','linewidth',linewidth)

set(gca,'linewidth',linewidth,'fontsize',fontsize,'Fontweight','bold')

ylabel('Reappraisal Success')
h = xlabel('Map Activation (Beta value)','FontSize',fontsize,'Fontweight','bold');
get(h)
h = ylabel('Reappraisal Success','FontSize',fontsize,'Fontweight','bold');
get(h)

title('Common Appraisal')
%%%%%%%%%%%%%%%%%
subplot(1,4,1)
scatter(mean(Reg_Neg.dat(indexReappraisalOnly,:),1)',Success,dotsize,dotcolor,'.');
h=lsline

set(h(1),'color','r','linewidth',linewidth)

set(gca,'linewidth',linewidth,'fontsize',fontsize,'Fontweight','bold')

ylabel('Reappraisal Success')
h = xlabel('Map Activation (Beta value)','FontSize',fontsize,'Fontweight','bold');
get(h)
h = ylabel('Reappraisal Success','FontSize',fontsize,'Fontweight','bold');
get(h)

title('Reappraisal Only')
%%%%%%%

subplot(1,4,3)
scatter(mean(Reg_Neg.dat(indexLookOnly,:),1)',Success,dotsize,dotcolor,'.');
h=lsline

set(h(1),'color','#FF594C','linewidth',linewidth)

set(gca,'linewidth',linewidth,'fontsize',fontsize,'Fontweight','bold')

ylabel('Reappraisal Success')
h = xlabel('Map Activation (Beta value)','FontSize',fontsize,'Fontweight','bold');
get(h)
h = ylabel('Reappraisal Success','FontSize',fontsize,'Fontweight','bold');
get(h)

title('Un-Modifiable Emotion generation')

%%%%%%%%%%%%%%%%
subplot(1,4,4)

scatter(mean(Reg_Neg.dat(indexReappraisal_D,:),1)',Success,dotsize,dotcolor,'.');
h=lsline

set(h(1),'color','#FF594C','linewidth',linewidth)

set(gca,'linewidth',linewidth,'fontsize',fontsize,'Fontweight','bold')

ylabel('Reappraisal Success')
h = xlabel('Map Activation (Beta value)','FontSize',fontsize,'Fontweight','bold');
get(h)
h = ylabel('Reappraisal Success','FontSize',fontsize,'Fontweight','bold');
get(h)

title('Modifiable Emotion generation')


Reappraisal_D_pattern=Reg_Neg.dat(indexReappraisalOnly,:);

for i=1:size(Reg_rating,1)
    for j=1:size(Reg_rating,1)
        [R P]=corr(Reappraisal_D_pattern(:,i),Reappraisal_D_pattern(:,j));
        Rvalue(i,j)=R;
        d=Success(i)-Success(j);
        dvalue(i,j)=d;
        
    end
end
figure
imagesc(Rvalue)
set(gca,'clim',[0.1 0.2])

figure
imagesc(dvalue)

k=1;
for i=1:size(Reg_rating,1)
    for j=i+1:size(Reg_rating,1)
        Rvector(k)=Rvalue(i,j);
        dvector(k)=dvalue(i,j);
        k=k+1;
    end
end

[R P]=corr(Rvector',dvector','type','spearman');

figure
scatter(Rvector',dvector')

figure
scatter(rand(1,5),rand(1,5),500,'r','.')
hold on
scatter(rand(1,5),rand(1,5),500,'b','.')
legend('Study 1','Study 2')
set(gca,'fontsize',12,'fontweight','bold')


montage
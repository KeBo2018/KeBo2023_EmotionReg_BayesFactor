clear
load('C:\Users\KeBo\Documents\GitHub\KeBo2023_EmotionReg_BayesFactor\Final_SystemComponentMap_Consensus_AfterClusterControl.mat')

load('C:\Users\KeBo\Dropbox (Dartmouth College)\2021_Ke_Bo_reappraisal_Gianaros\Data\fMRI_data\AHAB_FullData_Meta.mat')

RegNeg_AHAB=image_math(Whole_Reg,Whole_Neg,'minus')
load('C:\Users\KeBo\Dropbox (Dartmouth College)\2021_Ke_Bo_reappraisal_Gianaros\Data\fMRI_data\PIP_FullData_Meta.mat')

RegNeg_PIP=image_math(Whole_Reg,Whole_Neg,'minus')

Reg_rating_AHAB=table2array(RegNeg_AHAB.metadata_table(:,9))
Neg_rating_AHAB=table2array(RegNeg_AHAB.metadata_table(:,10))
Neu_rating_AHAB=table2array(RegNeg_AHAB.metadata_table(:,11))

Reg_rating_PIP=table2array(RegNeg_PIP.metadata_table(:,9))
Neg_rating_PIP=table2array(RegNeg_PIP.metadata_table(:,10))
Neu_rating_PIP=table2array(RegNeg_PIP.metadata_table(:,11))


Success_AHAB=Neg_rating_AHAB-Reg_rating_AHAB; 
EmoAct_AHAB=Neg_rating_AHAB-Neu_rating_AHAB;


EmoAct_PIP=Neg_rating_PIP-Neu_rating_PIP;
Success_PIP=Neg_rating_PIP-Reg_rating_PIP;

Success=[Success_AHAB;Success_PIP];
Reg_Neg=image_math(RegNeg_AHAB,RegNeg_PIP,'concatenate')


%%%%%%%%%%Univariate%%%%%
[R(1) P(1)]=corr(mean(Reg_Neg.dat(indexOverlap,:),1)',Success);
[R(2) P(2)]=corr(mean(Reg_Neg.dat(indexReappraisalOnly,:),1)',Success);
[R(3) P(3)]=corr(mean(Reg_Neg.dat(indexReappraisal_D,:),1)',Success);
[R(4) P(4)]=corr(mean(Reg_Neg.dat(indexLookOnly,:),1)',Success);

%%%%%%%

figure
dotcolor1=[0.3010 0.7450 0.9330];
dotcolor2=[124 141 204]/256;
dotsize=500;
linewidth=3;
fontsize=12;
subplot(1,4,2)
scatter(mean(Reg_Neg.dat(indexOverlap,:),1)',Success,dotsize,dotcolor1,'.');
h=lsline
hold on
scatter(mean(Reg_Neg.dat(indexOverlap,183:end),1)',Success(183:end),dotsize,dotcolor2,'.');

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
hold on
scatter(mean(Reg_Neg.dat(indexReappraisalOnly,183:end),1)',Success(183:end),dotsize,dotcolor2,'.');

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
hold on
scatter(mean(Reg_Neg.dat(indexLookOnly,183:end),1)',Success(183:end),dotsize,dotcolor2,'.');

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
hold on
scatter(mean(Reg_Neg.dat(indexReappraisal_D,183:end),1)',Success(183:end),dotsize,dotcolor2,'.');

set(h(1),'color','#FF594C','linewidth',linewidth)

set(gca,'linewidth',linewidth,'fontsize',fontsize,'Fontweight','bold')

ylabel('Reappraisal Success')
h = xlabel('Map Activation (Beta value)','FontSize',fontsize,'Fontweight','bold');
get(h)
h = ylabel('Reappraisal Success','FontSize',fontsize,'Fontweight','bold');
get(h)

title('Modifiable Emotion generation')


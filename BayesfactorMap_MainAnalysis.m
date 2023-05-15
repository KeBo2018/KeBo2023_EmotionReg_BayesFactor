clear
%%%% Select study 1 %%%
load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\AHAB_WithoutRest3_Whole_Beta.mat')

%%%%Select study 2 %%%
% load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\PIP_WithoutRest3_Whole_Beta.mat') 



RegMinusNeg = image_math(Whole_Reg, Whole_Neg,'minus' )  %%% Difference betweem Reappraisal and Negative look 
NegMinusNeu = image_math(Whole_Neg, Whole_Neu,'minus' )  %%% Negative vs Neutral

wholeSub=size(RegMinusNeg.dat,2);



tRegNeg=ttest(RegMinusNeg)      %%% T map generate

%%%% visualize T maps  %%%%
% t1 = threshold(tRegNeg, .05, 'fdr')
% orthviews(t1)
% RegNeg_FDR_P=t1.threshold;

% 
% create_figure('montage'); axis off;
% montage(t1,'compact');
% % drawnow, snapnow

tNegNeu= ttest(NegMinusNeu)
% t2 = threshold(tNegNeu, .05, 'fdr')
% 
% t2 = threshold(tNegNeu,[2.5 Inf], 'raw-between')
% orthviews(t2)
% create_figure('montage'); axis off;
% montage(t2,'compact');
% figure
% surface(t2)


%%%%%Compute Bayes factor
BF_RegNeg=estimateBayesFactor(tRegNeg,'t');   %%%%Bayes factor, note the output is 2*log(BF), BF=10 correspond to 2*log(bf)=4.6, BF=0.1 corresponding to -4,6
BF_NegNeu=estimateBayesFactor(tNegNeu,'t');


%  Visualize BF maps
% figure
% 
% t2 = threshold(BF_RegNeg, [-4.6 4.6], 'raw-outside')
% create_figure('montage'); axis off;
% montage(t2,'compact');
% figure
% surface(t2)
% 
% 
% t2 = threshold(BF_NegNeu, [-4.6 4.6], 'raw-outside')
% create_figure('montage'); axis off;
% montage(t2,'compact');
% figure
% surface(t2)
%%%%%%%%%%%%%%% Testing FDR threshold to BF threshold %%%%%%%%

%%%%test corresponding T for BF=10
sample_size=176;
scaling_factor=0.707;
myfun = @(T)t1smpbf(T,sample_size,scaling_factor);

% Define the target output value
BF_Threshold=10;
target = BF_Threshold;

% Use fzero to find the parameter that gives an output of 10
x0 = 1; % Initial guess for the parameter
x = fzero(@(x) myfun(x) - target, x0);

Threshold_T=x;

t=x;
p = 2 * (1 - tcdf(abs(t), sample_size));
%%%%%%%%%%%%%%%%%


Threshold1=4.6;   %%%% equal to bf=10
Threshold2=4.6;
tbase1=ttest(Whole_Neg);   %%%%% single condition activation for look negative
tbase2=ttest(Whole_Reg);    %%%%% single condition activation for look negative


%%%%%%%%%%%%%% Common appraisal %%%%%%
BF_tstat1=BF_RegNeg;
for i=1:242868
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) >Threshold2 && tRegNeg.dat(i)>0 && tNegNeu.dat(i)>0  && tbase1.dat(i)>0
       BF_tstat1.dat(i)=1;
   else
       BF_tstat1.dat(i)=0;
   end
       
end

Vcount_Overlap=length(find(BF_tstat1.dat==1));  %%%Voxel count of 'common appraisal '

%%%%%%%%%%% Reappraisal only %%%%%%%%%%%%%
BF_tstat2=BF_RegNeg;
for i=1:242868
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) <(-1*Threshold2) && tRegNeg.dat(i)>0 && tbase2.dat(i)>0
       BF_tstat2.dat(i)=1;
   else
       BF_tstat2.dat(i)=0;
   end
       
end

Vcount_ReappraisalOnly=length(find(BF_tstat2.dat==1))

%%%%%%% Non-modifiable emotion %%%%%
BF_tstat3=BF_RegNeg;
for i=1:242868
   if BF_RegNeg.dat(i) <-1*Threshold1 && BF_NegNeu.dat(i) >Threshold2 && tNegNeu.dat(i)>0 &&  tbase1.dat(i)>0
       BF_tstat3.dat(i)=1;
   else
       BF_tstat3.dat(i)=0;
   end
       
end

Vcount_LookOnly=length(find(BF_tstat3.dat==1))
%%%%%%%% Modifiable emotion %%%%%%%%%%%%


BF_tstat4=BF_RegNeg;
for i=1:242868
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) >(Threshold2) && tNegNeu.dat(i)>0 && tRegNeg.dat(i)<0 && tbase1.dat(i)>0
       BF_tstat4.dat(i)=1;
   else
       BF_tstat4.dat(i)=0;
   end
       
end

Vcount_Reappraisal_D_Only=length(find(BF_tstat4.dat==1))


%%%%%%%  Regions that are deactivated during emotion generation %%%%%%%%%
BF_tstat5=BF_RegNeg;
%%%%%%%%%%%%%%%
for i=1:242868
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) >Threshold2 && tRegNeg.dat(i)<0 && tNegNeu.dat(i)<0  
       BF_tstat5.dat(i)=1;
   else
       BF_tstat5.dat(i)=0;
   end
       
end

Vcount_CommonA_Decrease=length(find(BF_tstat5.dat==1))

BF_tstat6=BF_RegNeg;
for i=1:242868
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) <-1*Threshold2 && tRegNeg.dat(i)<0  
       BF_tstat6.dat(i)=1;
   else
       BF_tstat6.dat(i)=0;
   end
       
end

Vcount_ReappraisalOnly_Decrease=length(find(BF_tstat6.dat==1))

BF_tstat7=BF_RegNeg;
for i=1:242868
   if BF_RegNeg.dat(i) <-1*Threshold1 && BF_NegNeu.dat(i) >Threshold2 && tNegNeu.dat(i)<0  
       BF_tstat7.dat(i)=1;
   else
       BF_tstat7.dat(i)=0;
   end
       
end

Vcount_LookOnly_Decrease=length(find(BF_tstat7.dat==1))
% end
% fname='D:\CANlab_Note\Univariate analysis\Reappraisal_Pure_BF.nii'
% write(Reappraisal_Pure, 'fname', fname, 'overwrite');


indexReappraisalOnly=find(BF_tstat2.dat==1);
indexOverlap=find(BF_tstat1.dat==1);
indexLookOnly=find(BF_tstat3.dat==1);
indexReappraisal_D=find(BF_tstat4.dat==1)
indexCommonA_Decrease=find(BF_tstat5.dat==1)
indexReappraisalOnly_Decrease=find(BF_tstat6.dat==1)
indexLookOnly_Decrease=find(BF_tstat7.dat==1)

%%%% Average activation %%%%%%%

A(1)=mean(mean(Whole_Neu.dat(indexReappraisalOnly_Decrease,:)));
A(2)=mean(mean(Whole_Neg.dat(indexReappraisalOnly_Decrease,:)));
A(3)=mean(mean(Whole_Reg.dat(indexReappraisalOnly_Decrease,:)));

figure
x=[1;2;3]
bar(x(1),A(1),'g','LineWidth', 2)
hold on
bar(x(2),A(2),'b','LineWidth', 2)
hold on
bar(x(3),A(3),'r','LineWidth', 2)
set(gca,'linewidth',3,'Fontsize',20,'fontweight','bold')
ylabel('Activation')
% axis([0 4 0 1]) 
% xlabel({'Neutral','Negative look','Negative reappraisal'})


%%%%%%%%% compute within-participant standard deviation %%%%%%
S(1,:)=mean(Whole_Neu.dat(indexLookOnly,:),1)  %%%%%%%%% Mean activation across voxel in each condition %%%%
S(2,:)=mean(Whole_Neg.dat(indexLookOnly,:),1)
S(3,:)=mean(Whole_Reg.dat(indexLookOnly,:),1)

Grandmean=mean(mean(S,1)) %%%% Grand mean across condition

WithinS(1,:)=S(1,:)-mean(S,1)+Grandmean*ones(1,size(S,2));  %%% For each subject, substract their own mean across conditions and add up grand mean
WithinS(2,:)=S(2,:)-mean(S,1)+Grandmean*ones(1,size(S,2));  %%% Thus, subject-level variation is elimiated and within-subject variation is maintained
WithinS(3,:)=S(3,:)-mean(S,1)+Grandmean*ones(1,size(S,2));
 
errlow(1)=std(WithinS(1,:),1)/sqrt(size(S,2));
errlow(2)=std(WithinS(2,:),1)/sqrt(size(S,2));
errlow(3)=std(WithinS(3,:),1)/sqrt(size(S,2));
hold on
errhigh=errlow;
er = errorbar(x,A,errlow,errhigh,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

set(gca,'xtick',[])
names=['Neu';'Neg';'Reg']
set(gca,'xtick',[1:3],'xticklabel',names)

%%%%%%%%% Visualize system component maps %%%%%%

create_figure('montage'); axis off;
montage(BF_tstat2,'compact2');
drawnow, snapnow



%%%%%  Report regions in table %%%%%
cl = region(BF_tstat3);

% threshold based on extent of 8 vox or greater
num_vox_per_cluster = cat(1, cl.numVox);
cl(num_vox_per_cluster < 10) = [];

% Print a table, and return clusters with names attached in the cl structure
% This will separate positive and negative-valued voxels in each region

[clpos, clneg] = table(cl);
new_cl_with_labels = [clpos clneg];

% threshold based on extent of 50 vox or greater - just to have a set to
% display nicely:
num_vox_per_cluster = cat(1, new_cl_with_labels.numVox);
new_cl_with_labels(num_vox_per_cluster < 10) = [];

[r_pain, region_table_pain, table_legend_text_pain] = autolabel_regions_using_atlas(new_cl_with_labels); 
table(r_pain)

% Montage of the clusters, showing each
montage(new_cl_with_labels, 'colormap', 'regioncenters');
snapnow

montage(new_cl_with_labels,'colormap');
snapnow



%%%% Write map into disk %%%%%%


fname='D:\CANlab_Note\Univariate analysis\New\Bayes_LookOnly.nii'
write(BF_tstat3, 'fname', fname, 'overwrite');


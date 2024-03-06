clear
%%%% Select study 1, load beta maps %%%
% load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\AHAB_WithoutRest3_Whole_Beta.mat')

%%%% Select study 2, load beta maps %%%
load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\PIP_WithoutRest3_Whole_Beta.mat') 



RegMinusNeg = image_math(Whole_Reg, Whole_Neg,'minus' )  %%% Compute Regulate negative vs Look negative 
NegMinusNeu = image_math(Whole_Neg, Whole_Neu,'minus' )  %%% Compute Look negative vs Look neutral

%%% Generate T maps for Regulate negative vs Look negative 
tRegNeg=ttest(RegMinusNeg);  

%%% visualization  %%%
Example = threshold(tRegNeg, .05, 'fdr')
orthviews(Example)
create_figure('montage'); axis off;
montage(Example,'compact');
drawnow, snapnow
figure
surface(Example)

%%% Generate T maps for Look negative vs Look neutral 
tNegNeu= ttest(NegMinusNeu);

%%% visualization  %%%
Example = threshold(tNegNeu, .05, 'fdr')
orthviews(Example)
create_figure('montage'); axis off;
montage(Example,'compact');
figure
surface(Example)


%%%%% Compute Bayes factor from T value %%%%%%
BF_RegNeg=estimateBayesFactor(tRegNeg,'t');   %%%% logged Bayes factor, note the output is 2*log(BF), BF=10 correspond to 2*log(bf)=4.6, BF=0.1 corresponding to -4.6
BF_NegNeu=estimateBayesFactor(tNegNeu,'t');


%%Visualize BF maps
        figure
        
        Example = threshold(BF_RegNeg, [-4.6 4.6], 'raw-outside')
        create_figure('montage'); axis off;
        montage(Example,'compact');
        figure
        surface(Example)
        
        
        Example = threshold(BF_NegNeu, [-4.6 4.6], 'raw-outside')
        create_figure('montage'); axis off;
        montage(Example,'compact');
        figure
        surface(Example)

%%%%%%%%%%%%%%% Testing the corresponding p value to BF threshold %%%%%%%%

    %%%%test corresponding T for BF=10
    sample_size=size(NegMinusNeu.dat,2);
    scaling_factor=0.707;
    myfun = @(T)t1smpbf(T,sample_size,scaling_factor);
    
    % Define the target output value
    BF_Threshold=10; %Set it to 0.1 if you want to test the lower bound
    target = BF_Threshold;
    
    % Use fzero to find the parameter that gives an output of 10
    x0 = 1; % Initial guess for the parameter
    x = fzero(@(x) myfun(x) - target, x0);
    
    Threshold_T=x;
    
    t=x;
    p = 2 * (1 - tcdf(abs(t), sample_size));
%%%%%%%%%%%%%%%%%

%%%%%%% System identification approach starts here %%%%%

Threshold1=4.6;   %%%% In Reg vs Neg contrast, set the threshold for 2log(bf) that corresponds to bf=10
Threshold2=4.6;   %%%% In Neg vs Neu contrast, set the threshold for 2log(bf) that corresponds to bf=10

tbase1=ttest(Whole_Neg);   %%%%% single condition activation for look negative
tbase2=ttest(Whole_Reg);    %%%%% single condition activation for look negative


%%%%%%%%%%%%%% Define Common appraisal, make the map binary %%%%%%
BF_tstat1=BF_RegNeg;  %%Map initialization. The values will be changed in subsequent analysis

%Search across the whole brain, ‘Reappraise only’ voxels have BF > 10 and t > 0 for ‘Reappraisal’ and BF < 1/10 for ‘Emotion generation’. 
%In addition, these voxels have positive activation in the ' Look negative' condition.
for i=1:size(BF_tstat1.dat,1)
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) >Threshold2 && tRegNeg.dat(i)>0 && tNegNeu.dat(i)>0  && tbase1.dat(i)>0
       BF_tstat1.dat(i)=1;  %% Identified voxles are labeled as one
   else
       BF_tstat1.dat(i)=0;  %% Other voxels are labeled as zero
   end
       
end

Vcount_Overlap=length(find(BF_tstat1.dat==1));  %%%Voxel count of 'common appraisal '

%%%%%%%%%%% Reappraisal only %%%%%%%%%%%%%
BF_tstat2=BF_RegNeg;  %%Map initialization. The values will be changed in subsequent analysis

%‘Reappraise only’ voxels have BF > 10 and t>0 for Reappraisal and BF < 1/10 for Emotion Generation. 
% And these voxels have positive activation in ‘Regulate negative’ condition.

for i=1:size(BF_tstat2.dat,1)
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) <(-1*Threshold2) && tRegNeg.dat(i)>0 && tbase2.dat(i)>0
       BF_tstat2.dat(i)=1;
   else
       BF_tstat2.dat(i)=0;
   end
       
end

Vcount_ReappraisalOnly=length(find(BF_tstat2.dat==1))

%%%%%%% Non-modifiable emotion %%%%%
BF_tstat3=BF_RegNeg;
%%'Non-modifiable emotion-generation' voxels have BF > 10 and t>0 for Emotion Generation and BF < 1/10 for Reappraisal. 
%%And these voxels have positive activation in %Look negative% condition.

for i=1:size(BF_tstat3.dat,1)
   if BF_RegNeg.dat(i) <-1*Threshold1 && BF_NegNeu.dat(i) >Threshold2 && tNegNeu.dat(i)>0 &&  tbase1.dat(i)>0
       BF_tstat3.dat(i)=1;
   else
       BF_tstat3.dat(i)=0;
   end
       
end

Vcount_LookOnly=length(find(BF_tstat3.dat==1))
%%%%%%%% Modifiable emotion %%%%%%%%%%%%
BF_tstat4=BF_RegNeg;
%%'Modifiable emotion-generation% voxels have BF > 10 and t>0 for Emotion Generation and BF > 10 for Reappraisal with a negative effect (t<0, reduced activity during reappraisal). 
%%And these voxels have positive activation in %Look negative% condition.


for i=1:size(BF_tstat4.dat,1)
   if BF_RegNeg.dat(i) >Threshold1 && BF_NegNeu.dat(i) >(Threshold2) && tNegNeu.dat(i)>0 && tRegNeg.dat(i)<0 && tbase1.dat(i)>0
       BF_tstat4.dat(i)=1;
   else
       BF_tstat4.dat(i)=0;
   end
       
end

Vcount_Reappraisal_D=length(find(BF_tstat4.dat==1))



%%%% Save the indexes of the BF maps %%%
indexReappraisalOnly=find(BF_tstat2.dat==1);
indexOverlap=find(BF_tstat1.dat==1);
indexLookOnly=find(BF_tstat3.dat==1);
indexReappraisal_D=find(BF_tstat4.dat==1)

%%%% Save the BF maps using official names%%%
CommonAppraisal=BF_tstat1;
ReappraisalOnly=BF_tstat2;
NonModifiableEmo=BF_tstat3;
ModifiableEmo=BF_tstat4;

%%%% Average activation, customize index to view different maps %%%%%%%

A(1)=mean(mean(Whole_Neu.dat(indexReappraisalOnly,:)));
A(2)=mean(mean(Whole_Neg.dat(indexReappraisalOnly,:)));
A(3)=mean(mean(Whole_Reg.dat(indexReappraisalOnly,:)));

figure
x=[1;2;3]
bar(x(1),A(1),'g','LineWidth', 2)
hold on
bar(x(2),A(2),'b','LineWidth', 2)
hold on
bar(x(3),A(3),'r','LineWidth', 2)
set(gca,'linewidth',3,'Fontsize',20,'fontweight','bold')
ylabel('Activation')


%%%%%%%%% compute within-participant standard deviation, customize index to view different maps %%%%%%
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

%%%%%%%%% Visualize system component maps, customize map names to view different maps %%%%%%

create_figure('montage'); axis off;
montage(BF_tstat2,'compact2');
drawnow, snapnow



%%%%%  Report regions in table,customize map names to view report from different maps. Note this result won't have cluster control %%%%%
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

% 
% fname='D:\CANlab_Note\Univariate analysis\New\Bayes_LookOnly.nii'
% write(BF_tstat3, 'fname', fname, 'overwrite');


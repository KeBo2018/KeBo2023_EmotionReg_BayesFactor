clear
%%%%%%%%% Cluster control %%%%%%%%

%% Load the BF maps, choose from study1, 2 and concensus map %%
% load('C:\Users\KeBo\Documents\GitHub\KeBo2023_EmotionReg_BayesFactor\SystemComponents_Study1_BeforeClusterControl.mat')  %%Load BF maps
% load('C:\Users\KeBo\Documents\GitHub\KeBo2023_EmotionReg_BayesFactor\SystemComponents_Study2_BeforeClusterControl.mat')  %%Load BF maps
load('C:\Users\KeBo\Documents\GitHub\KeBo2023_EmotionReg_BayesFactor\Final_SystemComponentMap_Consensus_BeforeClusterControl.mat')  %%Load BF maps

%%Save them in a easier name ^^ %%%%
BF_tstat1=CommonAppraisal;
BF_tstat2=ReappraisalOnly;
BF_tstat3=NonModifiableEmo;
BF_tstat4=ModifiableEmo;

Subregion1=region(BF_tstat1);
Subregion2=region(BF_tstat2);
Subregion3=region(BF_tstat3);
Subregion4=region(BF_tstat4);
k=1;
for i=1:size(Subregion1,2)  % Search across all clusters and only select the one that is larger than 15 voxels
    
    if (length(Subregion1(i).val))>14
        Subregion1_CC(k)=Subregion1(i);
        k=k+1;
    end
end
k=1;
for i=1:size(Subregion2,2)
    
    if (length(Subregion2(i).val))>14
        Subregion2_CC(k)=Subregion2(i);
        k=k+1;
    end
end
k=1;
for i=1:size(Subregion3,2)
    
    if (length(Subregion3(i).val))>14
        Subregion3_CC(k)=Subregion3(i);
        k=k+1;
    end
end
k=1;
for i=1:size(Subregion4,2)
    
    if (length(Subregion4(i).val))>14
        Subregion4_CC(k)=Subregion4(i);
        k=k+1;
    end
end





%% visualize the maps that are controled %%


BF_tstat1_CC=region2fmri_data(Subregion1_CC(1:end),BF_tstat1)
orthviews(BF_tstat1_CC)

figure
BF_tstat2_CC=region2fmri_data(Subregion2_CC(1:end),BF_tstat2)
orthviews(BF_tstat2_CC)

BF_tstat3_CC=region2fmri_data(Subregion3_CC(1:end),BF_tstat3)
orthviews(BF_tstat3_CC)

BF_tstat4_CC=region2fmri_data(Subregion4_CC(1:end),BF_tstat4)
orthviews(BF_tstat4_CC)

index1=find(BF_tstat1_CC.dat>0);
BF_tstat1_CC.dat(index1)=1;

index2=find(BF_tstat2_CC.dat>0);
BF_tstat2_CC.dat(index2)=1;

index3=find(BF_tstat3_CC.dat>0);
BF_tstat3_CC.dat(index3)=1;

index4=find(BF_tstat4_CC.dat>0);
BF_tstat4_CC.dat(index4)=1;

figure
orthviews(BF_tstat1_CC)

%% Save the index
indexOverlap=index1;
indexReappraisalOnly=index2;
indexLookOnly=index3;
indexReappraisal_D=index4;
CommonAppraisal_CC=BF_tstat1_CC;
ReappraisalOnly_CC=BF_tstat2_CC;
NonModifiableEmo_CC=BF_tstat3_CC;
ModifiableEmo_CC=BF_tstat4_CC;

figure
surface(ReappraisalOnly_CC)
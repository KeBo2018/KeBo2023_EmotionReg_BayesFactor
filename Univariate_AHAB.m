clear
load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\Gianaros_Train_data.mat')

load('D:\CANlab_Working\CANlab_Note\Emotion Regulation\drive-download-20210518T210113Z-001\Gianaros_Test_data.mat')

Whole_Reg=image_math(train_reg, test_reg,'concatenate')  %%%Concatenate training and testing data
Whole_Neg=image_math(train_neg, test_neg,'concatenate')
Whole_Neu=image_math(train_neu, test_neu,'concatenate')

RegMinusNeg = image_math(Whole_Reg, Whole_Neg,'minus' )  %%%Difference betweem Reappraisal and Negative look 
NegMinusNeu = image_math(Whole_Neg, Whole_Neu,'minus' )  %%%Negative vs Neutral

tRegNeg=ttest(RegMinusNeg)      %%% T map generate
t1 = threshold(tRegNeg, .05, 'fdr')
orthviews(t1)
RegNeg_FDR_P=t1.threshold;

tbase=ttest(Whole_Neu)
% create_figure('montage'); axis off;
% montage(t2);
% drawnow, snapnow

tNegNeu= ttest(NegMinusNeu)
t2 = threshold(tNegNeu, .05, 'fdr')
orthviews(t2)

NegNeu_FDR_P=t2.threshold;
BF_RegNeg=estimateBayesFactor(tRegNeg,'t');   %%%%Bayes factor

% orthviews(BF_RegNeg_th)


BF_NegNeu=estimateBayesFactor(tNegNeu,'t');
% BF_NegNeu_th=threshold(BF_NegNeu,[-1000 3], 'raw-outside')



BF_RegNeg_th=threshold(BF_RegNeg,[-1000 -2.2], 'raw-between')
BF_NegNeu_th=threshold(BF_NegNeu,[-1000 2.2], 'raw-outside')   %%%% Setting threshold for Bayes factor
% orthviews(BF_NegNeu_th)
% orthviews(BF_RegNeg_th)

Reappraisal_Pure=conjunction(BF_RegNeg_th, BF_NegNeu_th); %%%Finding conjunction regions

for BF_Threshold=6:6
Threshold=2*log(BF_Threshold);


BF_tstat1=BF_RegNeg;
for i=1:328798
   if BF_RegNeg.dat(i) >Threshold && BF_NegNeu.dat(i) >Threshold && tRegNeg.dat(i)>0 && tNegNeu.dat(i)>0  && tbase.dat(i)>0
       BF_tstat1.dat(i)=1;
   else
       BF_tstat1.dat(i)=0;
   end
       
end

Vcount_Overlap(BF_Threshold)=length(find(BF_tstat1.dat==1))
% 
% create_figure('montage'); axis off;
% montage(BF_tstat3);
% drawnow, snapnow
%%%%%%%%%%%%%%%%%%%%%%%%
BF_tstat2=BF_RegNeg;
for i=1:328798
   if BF_RegNeg.dat(i) >Threshold && BF_NegNeu.dat(i) <(-1*Threshold) && tRegNeg.dat(i)>0 && tbase.dat(i)>0
       BF_tstat2.dat(i)=1;
   else
       BF_tstat2.dat(i)=0;
   end
       
end

Vcount_ReappraisalOnly(BF_Threshold)=length(find(BF_tstat2.dat==1))

%%%%%%%%%%%%
BF_tstat3=BF_RegNeg;
for i=1:328798
   if BF_RegNeg.dat(i) <-1*Threshold && BF_NegNeu.dat(i) >Threshold && tNegNeu.dat(i)>0 &&  tbase.dat(i)>0
       BF_tstat3.dat(i)=1;
   else
       BF_tstat3.dat(i)=0;
   end
       
end

Vcount_LookOnly(BF_Threshold)=length(find(BF_tstat3.dat==1))

end
% fname='D:\CANlab_Note\Univariate analysis\Reappraisal_Pure_BF.nii'
% write(Reappraisal_Pure, 'fname', fname, 'overwrite');


indexReappraisalOnly=find(BF_tstat2.dat==1);
indexOverlap=find(BF_tstat1.dat==1);
indexLookOnly=find(BF_tstat3.dat==1);

A(1)=mean(mean(Whole_Neu.dat(indexOverlap,:)));
A(2)=mean(mean(Whole_Neg.dat(indexOverlap,:)));
A(3)=mean(mean(Whole_Reg.dat(indexOverlap,:)));

figure
x=[1;2;3]
bar(x(1),A(1),'g','LineWidth', 2)
hold on
bar(x(2),A(2),'b','LineWidth', 2)
hold on
bar(x(3),A(3),'r','LineWidth', 2)
set(gca,'linewidth',3,'Fontsize',20,'fontweight','bold')
ylabel('Activation')
axis([0 4 0 0.5])
% xlabel({'Neutral','Negative look','Negative reappraisal'})

%  xaxis('Neutral')
%%%%%%%%%%%%% T map %%%
 tRegNeg=ttest(RegMinusNeg)      %%% T map generate
t1 = threshold(tRegNeg, .05, 'fdr')
orthviews(t1)

% create_figure('montage'); axis off;
% montage(t2);
% drawnow, snapnow

tNegNeu= ttest(NegMinusNeu)
t2 = threshold(tNegNeu,[0 1], 'raw-between')
orthviews(BF_tstat1)

Reappraisal_Pure=conjunction(t1, t2,1);


create_figure('montage'); axis off;
montage(BF_tstat3);
drawnow, snapnow


%%%%%%%%%%%%  encoding decoding model
train_neg_Overlap_Mat=train_neg.dat(indexOverlap,:);
test_neg_Overlap_Mat=test_neg.dat(indexOverlap,:);
train_reg_Overlap_Mat=train_reg.dat(indexOverlap,:);
test_reg_Overlap_Mat=test_reg.dat(indexOverlap,:);
train_neu_Overlap_Mat=train_neu.dat(indexOverlap,:);
test_neu_Overlap_Mat=test_neu.dat(indexOverlap,:);


train_neg_LookOnly_Mat=train_neg.dat(indexLookOnly,:);
test_neg_LookOnly_Mat=test_neg.dat(indexLookOnly,:);
train_reg_LookOnly_Mat=train_reg.dat(indexLookOnly,:);
test_reg_LookOnly_Mat=test_reg.dat(indexLookOnly,:);
train_neu_LookOnly_Mat=train_neu.dat(indexLookOnly,:);
test_neu_LookOnly_Mat=test_neu.dat(indexLookOnly,:);



train_neg_ReappraisalOnly_Mat=train_neg.dat(indexReappraisalOnly,:);
test_neg_ReappraisalOnly_Mat=test_neg.dat(indexReappraisalOnly,:);
train_reg_ReappraisalOnly_Mat=train_reg.dat(indexReappraisalOnly,:);
test_reg_ReappraisalOnly_Mat=test_reg.dat(indexReappraisalOnly,:);
train_neu_ReappraisalOnly_Mat=train_neu.dat(indexReappraisalOnly,:);
test_neu_ReappraisalOnly_Mat=test_neu.dat(indexReappraisalOnly,:);

%%%%%%%%%%


fname='D:\CANlab_Note\Univariate analysis\New\Bayes_LookOnly.nii'
write(BF_tstat3, 'fname', fname, 'overwrite');


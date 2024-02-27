# KeBo2023_EmotionReg_BayesFactor

Hi, Welcome to my Github page for code back up related to the publication 'A systems-identification approach using Bayes factors to deconstruct the brain bases of emotion regulation' in Nature Neuroscience.

Here are the brief workflow for the project:

Notice: Most of the code is based on Canlabtools. Please download canlabtools from https://github.com/canlab, and add Canlabcore and Neuroimaging_Pattern_Masks in your Matlab path.

1. The first level beta maps are stored in Neurovault: https://identifiers.org/neurovault.collection:16266. We also have ojbect level data (fMRI_data), which is the basic object used in Canlabtool. For the analysis in the paper, we need two types of beta contrast maps: Look negative vs Look neutral; and Regulate vs Look negative. Once these images from all participants are loaded, we can use the **'BayesfactorMap_MainAnalysis.m'** to create Bayes factor maps and perform corresponding system identificantion analysis for each individual datasests. The system component maps from both datasets are combined using **CreateConsensusMap.m**. The cluster size is controled above 15 voxels in any system component maps using code 'Cluster_Control_SC'. After cluster control, one should expect to get the final system component maps as shown in the paper.

2. Figure 2,3 5: The system components maps are plotted on brain using Canlab visualization. See specific code in **BrainFigureInManucript_Figure2_3_5.mlx**.

3. Figure 4:Amygdala and PAG ROI analysis are done using **Figure4_AmygdalaActivation**

4. Figure 6:The system componets maps are compared with Neurosynth topic maps and visualized on a two-dimensional PCA axis. See code in **BF_Map_Neurosynth_Topic_Ke.mlx**

5. Figure 7:The system componets maps are compared with Neurotransmitter PET map. See code in **Neurotransimitter_Association.mlx**

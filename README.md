# tsa_icc

tsaCalc.m
Code for calculating Time Shift Analysis (TSA), cross-correlations between all voxels' timeseries and a reference timeseries.
Here, the average signal from the venus sinuses is used as the reference signal as in the included references, there is also an option to include the global average signal
Apart from a 3D mask of the BOLD signal lead or lag for each voxel, the code also calculates the percentage of lead or lag voxels in each of the 90 AAL cortical regions. 

iccCalc.m
Code for calculating Intrinsic Connectivity Contrast (ICC), global connectivity of any given given voxel with all others in the cortex. 
Apart from a 3D mask of the strenght of global connectivity for each voxel, two (high and low) binary masks are crated including voxels that are in the e.g. 15th and 85th percentiles

For details please check our papers: 

Kavroulakis, E., Simos, N.J., Maris, T.G., Zaganas, I., Panagiotakis, S., & Papadaki, E.
Evidence of Age-Related Hemodynamic and Functional Connectivity Impairment: A Resting State fMRI Study.
Frontiers in Neurology (2021).

Antypa, D., Simos, N.J., Kavroulakis, E., Bertsias, G., Fanouriakis, A., Sidiropoulos, P., Boumpas, D., & Papadaki, E. 
Anxiety and depression severity in neuropsychiatric SLE are associated with perfusion and functional connectivity changes of the frontolimbic neural circuit: A resting-state f(unctional) MRI study. 
Lupus Science and Medicine (2021).

Papadaki, E., Simos, N.J., Kavroulakis, E. et al.
Converging evidence of impaired brain function in systemic lupus erythematosus: changes in perfusion dynamics and intrinsic functional connectivity.
Neuroradiology (2022).


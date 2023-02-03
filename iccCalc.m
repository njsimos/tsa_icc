% Code for calculating Intrinsic Connectivity Contrast (ICC), global connectivity of any given given voxel with all others in the cortex. 
% Apart from a 3D mask of the strenght of global connectivity for each voxel, two (high and low) binary masks are crated including voxels that are in the e.g. 15th and 85th percentiles
% For details please check our papers: 

% Kavroulakis, E., Simos, N.J., Maris, T.G., Zaganas, I., Panagiotakis, S., & Papadaki, E.
% Evidence of Age-Related Hemodynamic and Functional Connectivity Impairment: A Resting State fMRI Study.
% Frontiers in Neurology (2021).

% Antypa, D., Simos, N.J., Kavroulakis, E., Bertsias, G., Fanouriakis, A., Sidiropoulos, P., Boumpas, D., & Papadaki, E. 
% Anxiety and depression severity in neuropsychiatric SLE are associated with perfusion and functional connectivity changes of the frontolimbic neural circuit: A resting-state f(unctional) MRI study. 
% Lupus Science and Medicine (2021).

% Papadaki, E., Simos, N.J., Kavroulakis, E. et al.
% Converging evidence of impaired brain function in systemic lupus erythematosus: changes in perfusion dynamics and intrinsic functional connectivity.
% Neuroradiology (2022).

% Author: Nicholas John Simos, 2021

clear all
close all
clc

%%
% User Declarations, mandatory
dataDir = ''; % Direcotry of 4D functional data
dataName = ''; % e.g. zc_bd_C_4D.nii
voxel_size = [2 2 2];
origin = [40 57 26];
lowPrc = 15; % High and low percentile thresholds for creating binary high and low ICC mask
highPrc = 85; 
%%
datatype = 16;
% AAL entire mask for extracting cortex
aal3D_mask = load_nii('AAL_whole_centered_mask.nii').img;
numOfvox = sum(nonzeros(aal3D_mask));

sub4D_data = single(load_nii(strcat(dataDir, '\', dataName)).img);
numOfTp = size(sub4D_data, 4);
aalVoxTsSub = zeros(numOfvox, numOfTp);
aal3D_coord = zeros(numOfvox, 3);
voxCount = 1;
for i = 1:size(sub4D_data, 1)
    for j = 1:size(sub4D_data, 2)
        for k = 1:size(sub4D_data, 3)
            if aal3D_mask(i,j,k) == 1
                aalVoxTsSub(voxCount,:) = sub4D_data(i,j,k,:);
                aal3D_coord(voxCount,:) = [i,j,k];
                voxCount = voxCount + 1;
            end
        end    
    end
end

ICC = intrinsic_connectivity_contrast(aalVoxTsSub');
ICC_norm = normalize(ICC, 'range');
ICC_sqrt = sqrt(ICC);

reconstr3D_ICC_sqrt = zeros(size(sub4D_data, 1:3));

%reconstructing a 3D image from the 1D ICC vector (contains all the icc values for voxels inside AAL)
for voxCount = 1:length(aalVoxTsSub)
    reconstr3D_ICC_sqrt(aal3D_coord(voxCount,1),aal3D_coord(voxCount,2),aal3D_coord(voxCount,3),:) = ICC_sqrt(voxCount);
end

reconstr3D_ICC_sqrt = reconstr3D_ICC_sqrt.*double(aal3D_mask);

lowPrcIcc = prctile(reshape(reconstr3D_ICC_sqrt(aal3D_mask==1),1,[]), lowPrc);
highPrcIcc = prctile(reshape(reconstr3D_ICC_sqrt(aal3D_mask==1),1,[]), highPrc);

reconstr3D_ICC_loBin = zeros(size(reconstr3D_ICC_sqrt));
reconstr3D_ICC_loBin(reconstr3D_ICC_sqrt<=lowPrcIcc) = 1;
reconstr3D_ICC_loBin(reconstr3D_ICC_sqrt==0) = 0;
reconstr3D_ICC_loBin = uint8(reconstr3D_ICC_loBin.*double(aal3D_mask));

reconstr3D_ICC_hiBin = zeros(size(reconstr3D_ICC_sqrt));
reconstr3D_ICC_hiBin(reconstr3D_ICC_sqrt>=highPrcIcc) = 1;
reconstr3D_ICC_hiBin = uint8(reconstr3D_ICC_hiBin);

cd(dataDir)
%save subject icc maps
nii = make_nii(reconstr3D_ICC_sqrt, voxel_size, origin, datatype);
save_nii(nii, 'iccVals_sqrt.nii')

%save subject low and high icc maps
nii = make_nii(reconstr3D_ICC_loBin, voxel_size, origin, datatype);
save_nii(nii, 'binLowIcc.nii')
nii = make_nii(reconstr3D_ICC_hiBin, voxel_size, origin, datatype);
save_nii(nii, '_binHighIcc.nii')


% figure
% montage(reconstr3D_ICC_hiBin, 'DisplayRange', [-1 2])
% figure
% montage(reconstr3D_ICC_loBin, 'DisplayRange', [-1 2])


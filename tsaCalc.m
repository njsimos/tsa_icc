% Code for calculating Time Shift Analysis (TSA), cross-correlations between all voxels' timeseries and a reference timeseries.
% Here, the average signal from the venus sinuses is used as the reference signal as in the included references, there is also an option to include the global average signal
% Apart from a 3D mask of the BOLD signal lead or lag for each voxel, the code also calculates the percentage of lead or lag voxels in each of the 90 AAL cortical regions. 
% These are then stored as the two rows of the automatically saved file 'roiNegPosDelayPerc' (.mat and .xlsx), percentage of lag in first and lead in second rows
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
% User Declarations: Mandatory!
dataDir = ''; % Direcotry of 4D functional data
dataName = ''; % e.g. zc_bd_C_4D.nii
masksDir = ''; % Directory of AAL masks, inlcuded with this script
voxel_size = [2 2 2];
origin = [40 57 26];
TR = 2.32;
lagTRs = 3; % How many TRs to calculate lags for, e.g. [-3 3] range
%%

cd(dataDir)
datatype = 16; %denoting float32 otherwise 4-int16

% Venus sinus mask for extracting reference timeseries
venat3D_mask = load_nii('VENAT_PartialVolume_cent_thresh_2mm.nii').img;
numOfvoxVen = sum(nonzeros(venat3D_mask));

% AAL entire mask for extracting cortex
glob3D_mask = load_nii('AAL_whole_centered_mask.nii').img;
numOfvox = sum(nonzeros(glob3D_mask));

sub4D_data = single(load_nii(strcat(dataDir, '\', dataName)).img);
numOfTp = size(sub4D_data, 4);
aalVoxTsSub = zeros(numOfvox, numOfTp);
venVoxTsSub = zeros(numOfvoxVen, numOfTp);
aal3D_coord = zeros(numOfvox, 3);
voxCount = 1;
voxCountVen = 1;
for i = 1:size(sub4D_data, 1)
    for j = 1:size(sub4D_data, 2)
        for k = 1:size(sub4D_data, 3)
            if glob3D_mask(i,j,k) == 1
                aalVoxTsSub(voxCount,:) = sub4D_data(i,j,k,:);
                aal3D_coord(voxCount,:) = [i,j,k];
                voxCount = voxCount + 1;
            end
            if venat3D_mask(i,j,k) == 1
                venVoxTsSub(voxCountVen,:) = sub4D_data(i,j,k,:);
                voxCountVen = voxCountVen + 1;
            end

        end    
    end
end

r = zeros(1, lagTRs*2+1);
lags = zeros(1, lagTRs*2+1);

globalTsVen = mean(venVoxTsSub,1);
bestLagFitVen = zeros(1, numOfvox);
for tsNum = 1:numOfvox

    [r, lags] = xcorr(aalVoxTsSub(tsNum, :), globalTsVen, 3, 'coeff');
    [B, indices] = sort(r, 'descend');
    bestLagFitVen(tsNum) = lags(indices(1));

end

reconstr3D_lagVen = zeros(size(sub4D_data, 1:3));
for voxCount = 1:length(aalVoxTsSub)
    reconstr3D_lagVen(aal3D_coord(voxCount,1),aal3D_coord(voxCount,2),aal3D_coord(voxCount,3),:) = bestLagFitVen(voxCount)*TR;
end

% Binary positive and negative masks
tsa_pos_bin = zeros(size(reconstr3D_lagVen));
tsa_neg_bin = zeros(size(reconstr3D_lagVen));

tsa_pos_bin(reconstr3D_lagVen>=2*TR) = 1;
tsa_neg_bin(reconstr3D_lagVen<=-2*TR) = 1;

tsa_pos_bin = uint8(tsa_pos_bin).*glob3D_mask;
tsa_neg_bin = uint8(tsa_neg_bin).*glob3D_mask;

nii = make_nii(reconstr3D_lagVen, voxel_size, origin);
save_nii(nii, 'voxDelTimeFl_ven.nii')

nii = make_nii(tsa_pos_bin, voxel_size, origin);
save_nii(nii, 'tsa_pos_bin.nii')

nii = make_nii(tsa_neg_bin, voxel_size, origin);
save_nii(nii, 'tsa_neg_bin.nii')
    
    % Uncomment this section if TSA calculated with global time series delay is prefered
    %{
    globalTs = mean(aalVoxTsSub,1);
    bestLagFitGlob = zeros(1, numOfvox);
    for tsNum = 1:numOfvox
        [r, lags] = xcorr(aalVoxTsSub(tsNum, :), globalTs, 3, 'coeff');
        [B, indices] = sort(r, 'descend');
        bestLagFitGlob(tsNum) = lags(indices(1));
        
    end
    
    reconstr3D_lagGlob = zeros(size(sub4D_data, 1:3));
    for voxCount = 1:length(aalVoxTsSub)
        reconstr3D_lagGlob(aal3D_coord(voxCount,1),aal3D_coord(voxCount,2),aal3D_coord(voxCount,3),:) = bestLagFitGlob(voxCount)*TR;
    end
    
    nii = make_nii(reconstr3D_lagGlob, voxel_size, origin);
    save_nii(nii, 'voxDelTimeFl_glob.nii')
    %}
    
    % Display TSA values
    % montage(reconstr3D_lagVen, 'DisplayRange', [-6.96 6.96])
    
    
% TSA post-calculation metrics, ROI percentages 
roiNum = 1;
vec = ["2001","2002","2101","2102","2111","2112","2211","2212","2201","2202","2301","2302","2311","2312","2321","2322","2331","2332","2401","2402","2501","2502","2601","2602","2611","2612",...
"2701","2702","3001","3002","4001","4002","4011","4012","4021","4022","4101","4102","4111","4112","4201","4202","5001","5002","5011","5012","5021","5022","5101","5102","5201","5202","5301",...
"5302","5401","5402","6001","6002","6101","6102","6201","6202","6211","6212","6221","6222","6301","6302","6401","6402","7001","7002","7011","7012","7021","7022","7101","7102","8101","8102",...
"8111","8112","8121","8122","8201","8202","8211","8212","8301","8302"];
for roiName = vec
    AALmasks(roiNum).mask = cast(load_nii(strcat(masksDir, '\AAL1_rois_cent\AAL_', char(roiName), '.nii')).img, 'int16');
    roiNum = roiNum + 1;
end

load('AAL_roiNames.mat')

roiNegDelayPerc = zeros(1, length(AALmasks));
roiPosDelayPerc = zeros(1, length(AALmasks));
roiLength = zeros(1, length(AALmasks));

%masking the three different vox delay time matrices (raw, positive,
%negative) with the AAL ROIS
for roiNum = 1:length(AALmasks)
    voxCount = 1;
    voxCountPos = 1;
    voxCountNeg = 1;

    for i = 1:size(reconstr3D_lagVen, 1)
        for j = 1:size(reconstr3D_lagVen, 2)
            for k = 1:size(reconstr3D_lagVen, 3)
                if AALmasks(roiNum).mask(i,j,k) == 1
                    roiDelays(roiNum).val(voxCount) = reconstr3D_lagVen(i,j,k);

                    if reconstr3D_lagVen(i,j,k) > TR
                        roiPosDelays(roiNum).val(voxCount) = reconstr3D_lagVen(i,j,k);
                        voxCountPos = voxCountPos + 1;
                    else
                        roiPosDelays(roiNum).val(voxCount) = 0;
                    end
                    if reconstr3D_lagVen(i,j,k) < -TR
                        roiNegDelays(roiNum).val(voxCount) = -reconstr3D_lagVen(i,j,k);
                        voxCountNeg = voxCountNeg + 1;
                    else
                        roiNegDelays(roiNum).val(voxCount) = 0;
                    end

                    voxCount = voxCount + 1;
                end
            end    
        end
    end

    %calculation of number of positive and negative time delays for each ROI
    roiLength(roiNum) = length(roiDelays(roiNum).val);

    roiNegDelayPerc(roiNum) = 1 - (roiLength(roiNum)-voxCountNeg)/roiLength(roiNum);
    roiPosDelayPerc(roiNum) = 1 - (roiLength(roiNum)-voxCountPos)/roiLength(roiNum);

end

% Regions that displayed higher positive delay overall in their voxels
disp('Positive significant: ')
AAL_roiNames(roiPosDelayPerc>prctile(roiPosDelayPerc,80))

% Regions that displayed higher negative delay overall in their voxels
disp('Negative significant): ')
AAL_roiNames(roiNegDelayPerc>prctile(roiNegDelayPerc,80))

roiNegPosDelayPerc = cat(1, roiNegDelayPerc, roiPosDelayPerc);
cd(dataDir)
save('roiNegPosDelayPerc','roiNegPosDelayPerc')
xlswrite('roiNegPosDelayPerc',roiNegPosDelayPerc)



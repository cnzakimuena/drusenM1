
function call_drusenM1()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name - drusenM1
% Creation Date - 17th January 2021
% Author - Charles Belanger Nzakimuena
% Website - https://www.ibis-space.com/
%
% Description - 
%   The 'drusenM1' algorithm uses segmentation data to characterize the è
%   space between the retinal pigment epithelium (RPE) and Bruch’s 
%   membrane (BM).  Area and volume values corresponding to ETDRS subfieds 
%   are provided in table format.
%
% Example -
%		call_drusenM1()
%
% License - MIT
%
% Change History -
%                   17th January 2021 - Creation by Charles Belanger Nzakimuena
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('./subfunctions'))

%% list names of folders inside the patients folder

currentFolder = pwd;
patientsFolder = fullfile(currentFolder, 'processed');
myDir = dir(patientsFolder);
dirFlags = [myDir.isdir] & ~strcmp({myDir.name},'.') & ~strcmp({myDir.name},'..');
nameFolds = myDir(dirFlags);

%% for each 3x3 subfolder, turn segmented data into network graph

% get table row count
rowCount = 0;
for g = 1:numel(nameFolds)
    folder2 = fullfile(patientsFolder, nameFolds(g).name);
    patientDir2 = dir(fullfile(folder2, 'Results'));
    dirFlags2 = [patientDir2.isdir] & ~strcmp({patientDir2.name},'.') & ~strcmp({patientDir2.name},'..');
    subFolders2 = patientDir2(dirFlags2);
    rowCount = rowCount + numel(subFolders2);
end

col = zeros(rowCount,1);
colc = cell(rowCount,1);
% aTable = table(colc,col,col,col,col,col,col,col,col,col,col,...
%     'VariableNames',{'id' 'totalArea' 'region1' 'region2' 'region3' 'region4'...
%     'region5' 'jaccard' 'dice' 'sensitivity' 'specificity'});
aTable = table(colc,col,col,col,col,col,col,...
    'VariableNames',{'id' 'totalArea' 'region1' 'region2' 'region3' 'region4'...
    'region5'});
vTable = table(colc,col,col,col,col,col,col,...
    'VariableNames',{'id' 'totalVolume' 'region1' 'region2' 'region3' 'region4'...
    'region5'});

tableRow = 0;
for i = 1:numel(nameFolds)
    
    % assemble patient folder string
    folder = fullfile(patientsFolder, nameFolds(i).name);
    
    % add line to LOG
    disp(logit(folder, ['Initiating drusenM1; ' nameFolds(i).name ' folder']))
    
    patientDir = dir(fullfile(folder, 'Results'));
    dirFlags = [patientDir.isdir] & ~strcmp({patientDir.name},'.') & ~strcmp({patientDir.name},'..');
    subFolders = patientDir(dirFlags);

    for k = 1:numel(subFolders)
        
        nameFold = subFolders(k).name;
        scanType = nameFold(1:2);
        
        if strcmp(scanType, '3m')
            
            load(fullfile(folder,'Results', nameFold,'segmentation.mat'));
            load(fullfile(folder,'Results', nameFold,'scanInfo.mat'));   
            load(fullfile(folder,'Results', nameFold, 'ETDRS_grid','2DregionsETDRS.mat')); 
            load(fullfile(folder,'Results', nameFold, 'ETDRS_grid','3DregionsETDRS.mat')); 
            sizeRed = scanTag{2};
            
            % loop below prevents RPE bottom from falling below the BM
            for q = 1:size(RPEb,1)
                for qq = 1:size(RPEb,2)
                    if RPEb(q,qq) > lBM(q,qq)
                        RPEb(q,qq) = lBM(q,qq);
                    end
                end
            end
            dSpace = lBM-RPEb;

            % MODIFY THE MAP PROPORTIONS/ORIENTATION FOR EXPORT HERE
            numcols00 = round(size(dSpace, 1)*1536/300*sizeRed);
            numrows00 = round(size(dSpace, 2)*sizeRed);
            dSpace1 = imresize(dSpace,[numrows00 numcols00]);
                if strcmp(scanTag{1},'OD')
                    dSpace2 = flip(dSpace1,2);
                elseif strcmp(scanTag{1},'OS')
                    dSpace2 = flip(dSpace1,2);
                end
%           figure;imshow(dSpace2,[])      
            
            % x(:) transforms the array to a column vector
            dSpace_avg = mean(dSpace(:));
            dSpace_std = std(dSpace(:));
            
%             figure;imshow3D(volumeStruc,[],'plot',cat(3,lBM,RPEb,RVIf),'LineWidth',2)
%             bscanNum = 185;
%             chartColors1.c2 = rgb('RoyalBlue');
%             chartColors2.c2 = rgb('LightSkyBlue');
%             chartColors3.c2 = rgb('HotPink');
%             figure;imshow(volumeStruc(:,:,bscanNum),[])
%             hold on;
%             plot(lBM(bscanNum,:),'color',chartColors1.c2,'LineWidth',3);
%             plot(RPEb(bscanNum,:),'color',chartColors2.c2,'LineWidth',3);
% %             plot(RPEt(bscanNum,:),'color',chartColors1.c2,'LineWidth',2);
% %             plot(RVIf(bscanNum,:),'color',chartColors3.c2,'LineWidth',2);
            
            volumeMask = zeros(size(volumeStruc),'like',volumeStruc);
            
            % thresholds that detmines minimum dSpace depth that is kept as
            % the drusen mask, as a factor of 'dSpace_std'
            depthThresh = 0.5;
            absTresh = 20*sizeRed; 
            
            % iterating through each element of the dSpace 2D array, if one of the
            % elements exceeds the standard deviation, the space between lBM and RPEb
            % in volumeMask is assigned the value of 1
            for g = 1:size(dSpace,1)
                for gg = 1:size(dSpace,2)
                    if abs(dSpace(g,gg)) > (absTresh+depthThresh*dSpace_std)
                        if lBM(g,gg) > RPEb(g,gg)
                            volumeMask(RPEb(g,gg):lBM(g,gg),gg,g) = 1;
                            % % add 'elseif' statement to include sub-BM spaces in mask
                            % elseif RPEb(g,gg) > lBM(g,gg)
                            % volumeMask(lBM(g,ii):RPEb(g,ii),ii,g) = 1;
                        end
                    end
                end
            end
            % figure;imshow3D(volumeMask,[],'plot',cat(3,lBM,RPEb),'LineWidth',2)
%             bscanNum = 164;
%             chartColors1.c2 = rgb('RoyalBlue');
%             chartColors2.c2 = rgb('LightSkyBlue');
%             figure;imshow(volumeMask(:,:,bscanNum),[])
%             hold on;
%             plot(lBM(bscanNum,:),'color',chartColors1.c2,'LineWidth',2);
%             plot(RPEb(bscanNum,:),'color',chartColors2.c2,'LineWidth',2);
            
            % volumeMask orientation change to en-face direction
            enFace_Mask = [];
            for f = 1:size(volumeMask,1)
                enFace_im = mat2gray(reshape(volumeMask(f,:,:), [size(volumeMask, 2), size(volumeMask, 3)]));
                if strcmp(scanTag{1},'OS')
                    enFace_Mask(:,:,f) = imrotate(enFace_im,-90);
                elseif strcmp(scanTag{1},'OD')
                    enFace_Mask(:,:,f) = imrotate(enFace_im,-90);
                end
            end
            %figure;imshow3D(enFace_Mask,[])
            
            % check data orientation integrity
            maxMask = zeros(size(enFace_Mask, 1),size(enFace_Mask, 2));
            for h = 1:size(enFace_Mask, 1)
                for hh = 1:size(enFace_Mask, 2)
                    maxMask(h,hh) = max(enFace_Mask(h,hh,:));
                end
            end

%             % volumeStruc orientation change to en-face direction
%             enFace_Struc = [];
%             for ff = 1:size(volumeStruc,1)
%                 enFace_im2 = mat2gray(reshape(volumeStruc(ff,:,:), [size(volumeStruc, 2), size(volumeStruc, 3)]));
%                 if strcmp(scanTag{1},'OS')
%                     enFace_Struc(:,:,ff) = imrotate(enFace_im2,-90);
%                 elseif strcmp(scanTag{1},'OD')
%                     enFace_Struc(:,:,ff) = imrotate(enFace_im2,-90);
%                 end
%             end
%             %figure;imshow3D(enFace_Struc,[])
%             
%             % check data orientation integrity
%             maxStruc = zeros(size(enFace_Struc, 1),size(enFace_Struc, 2));
%             for d = 1:size(enFace_Struc, 1)
%                 for dd = 1:size(enFace_Struc, 2)
%                     maxStruc(d,dd) = max(enFace_Struc(d,dd,:));
%                 end
%             end
%             % figure;imshow(maxStruc,[])
 
            % modify 3D volume aspect ratio to true proportions
            % dimension decreased by reduction factor
            decNewDimA = round(size(enFace_Mask, 1)*1536/300*sizeRed);
            decNewDimB = round(size(enFace_Mask, 2)*sizeRed);
            % assuming same real distance in a-scan and volume direction (3mmx3mm)
            decNewDimV = round(size(enFace_Mask, 3)*sizeRed);
            % returns the volume B that has the number of rows, columns, and planes
            % specified by the three-element vector [numrows numcols numplanes].
            enFace_Mask_BW = imbinarize(imresize3(enFace_Mask,[decNewDimA decNewDimB decNewDimV]));
            %figure;imshow3D(enFace_Mask_BW,[])            
            
%             numcols = round(size(maxStruc, 1)*1536/300*sizeRed);
%             numrows = round(size(maxStruc, 2)*sizeRed);
%             maxStruc2 = imresize(maxStruc,[numrows numcols]);
%             figure;imshow(maxStruc2,[])      

            numcols2 = round(size(maxMask, 1)*1536/300*sizeRed);
            numrows2 = round(size(maxMask, 2)*sizeRed);
            maxMask_BW = imbinarize(imresize(maxMask,[numrows2 numcols2]));
            imwrite(maxMask_BW,fullfile([folder,'\Results\', nameFold, '\drusenMaskM1_BW' '.png']));
%             figure;imshow(maxMask_BW,[])           
% %             chartColors1.c2 = rgb('Black');
% %             chartColors2.c2 = rgb('LightSkyBlue');
% %             colormap([chartColors1.c2; chartColors2.c2]) % Blue = [0 0 1], White = [1 1 1]
            
%             % obtain validation parameters
%             if strcmp(patientsType,'AMD')
%                 valid_folder = 'C:\Users\...';
%                 valid_Mask = imread(fullfile([valid_folder,'\',nameFolds(i).name,'\',nameFold, '.png']));
%                 if size(valid_Mask,3)==3
%                     valid_Mask = rgb2gray(valid_Mask);
%                 end
%                 level = graythresh(valid_Mask);
%                 valid_MaskBW = imbinarize(valid_Mask,level);
%                 valid_MaskBW = imresize(valid_MaskBW,[numrows2 numcols2]);
%             elseif strcmp(patientsType,'normal') || strcmp(patientsType,'normal_U50') 
%                 valid_MaskBW = imbinarize(zeros(size(maxMask_BW)));
%             end
%             jaccard_coef = jaccard(valid_MaskBW,maxMask_BW);
%             dice_coef = dice(valid_MaskBW,maxMask_BW);
%             TP=0;FP=0;TN=0;FN=0;
%             for w=1:size(maxMask_BW, 1)
%                 for ww=1:400
%                     if(valid_MaskBW(w,ww)==1 && maxMask_BW(w,ww)==1)
%                         TP=TP+1;
%                     elseif(valid_MaskBW(w,ww)==0 && maxMask_BW(w,ww)==1)
%                         FP=FP+1;
%                     elseif(valid_MaskBW(w,ww)==0 && maxMask_BW(w,ww)==0)
%                         TN=TN+1;
%                     else
%                         FN=FN+1;
%                     end
%                 end
%             end
%             sensitivity = TP/(TP+FN);
%             specificity = TN/(TN+FP);
   
%             % data export (enFace_Mask_BW, structETDRS)
            dataExport2(folder, nameFold, enFace_Mask_BW, structETDRS)

            % 2D ETDRS regions setup
            disp('begin fundProfile_2D')
            aProfile = fundProfile_2D(maxMask_BW, regionsETDRS, sizeRed);
            disp('end fundProfile_2D')
            
%             aProfile = [aProfile jaccard_coef dice_coef sensitivity specificity];
            
            % 3D ETDRS regions setup
            disp('begin fundProfile_3D')
            vProfile = fundProfile_3D(enFace_Mask_BW, structETDRS, sizeRed);
            disp('end fundProfile_3D')

            % For left eye, ETDRS regions must be modified from OD nomenclature
            % to OS nomenclature
            if contains(nameFold, '_OS_')
                aRegion3 = aProfile(6);
                aRegion5 = aProfile(4);
                aProfile(4) = aRegion3;
                aProfile(6) = aRegion5;
                
                vRegion3 = vProfile(6);
                vRegion5 = vProfile(4);
                vProfile(4) = vRegion3;
                vProfile(6) = vRegion5;

            end
  
            tableRow = tableRow + 1;
            
            aTable{tableRow,'id'} = {nameFold};
            aTable{tableRow,'totalArea'} = aProfile(1);
            aTable{tableRow,'region1'}  = aProfile(2);
            aTable{tableRow,'region2'} = aProfile(3);
            aTable{tableRow,'region3'} = aProfile(4);
            aTable{tableRow,'region4'} = aProfile(5);
            aTable{tableRow,'region5'} = aProfile(6);
%             aTable{tableRow,'jaccard'} = aProfile(7);
%             aTable{tableRow,'dice'} = aProfile(8);
%             aTable{tableRow,'sensitivity'} = aProfile(9);
%             aTable{tableRow,'specificity'} = aProfile(10);
                        
            vTable{tableRow,'id'} = {nameFold};
            vTable{tableRow,'totalVolume'} = vProfile(1);
            vTable{tableRow,'region1'}  = vProfile(2);
            vTable{tableRow,'region2'} = vProfile(3);
            vTable{tableRow,'region3'} = vProfile(4);
            vTable{tableRow,'region4'} = vProfile(5);
            vTable{tableRow,'region5'} = vProfile(6);

        end
    end
    
end

fileName1 = fullfile(patientsFolder,'aTable.xls');
fileName2 = fullfile(patientsFolder,'vTable.xls');
writetable(aTable,fileName1)
writetable(vTable,fileName2)

disp(logit(patientsFolder,'Done drusenM1'))

            
            

function dataExport2(exportFolder, folderName, enFace_BW, sETDRS)

if ~exist(fullfile(exportFolder,'ExportFiles'), 'dir')
    mkdir(fullfile(exportFolder,'ExportFiles'));
end
if ~exist(fullfile(exportFolder,'ExportFiles', folderName), 'dir')
    mkdir(fullfile(exportFolder,'ExportFiles', folderName));
end

% 'enFace_BW' and 'sETDRS.regionsETDRS_3D{i}' volume data export
for i = 1:size(enFace_BW,3)
    imID = 'drusen_BW';
    im = enFace_BW(:,:,i);
    if ~exist(fullfile(exportFolder,'ExportFiles', folderName, imID), 'dir')
        mkdir(fullfile(exportFolder,'ExportFiles', folderName, imID));
    end
    imwrite(im, fullfile(exportFolder,'ExportFiles', folderName, imID, [imID '_' num2str(i, '%03.0f') '.png']))
end

for ii = 1:size(sETDRS.regionsETDRS_3D{1},3)
    imID = 'ETDRS_3D_1';
    im = sETDRS.regionsETDRS_3D{1}(:,:,ii);
    if ~exist(fullfile(exportFolder,'ExportFiles', folderName, imID), 'dir')
        mkdir(fullfile(exportFolder,'ExportFiles', folderName, imID));
    end
    imwrite(im, fullfile(exportFolder,'ExportFiles', folderName, imID, [imID '_' num2str(ii, '%03.0f') '.png']))
end
for k = 1:size(sETDRS.regionsETDRS_3D{2},3)
    imID = 'ETDRS_3D_2';
    im = sETDRS.regionsETDRS_3D{2}(:,:,k);
    if ~exist(fullfile(exportFolder,'ExportFiles', folderName, imID), 'dir')
        mkdir(fullfile(exportFolder,'ExportFiles', folderName, imID));
    end
    imwrite(im, fullfile(exportFolder,'ExportFiles', folderName, imID, [imID '_' num2str(k, '%03.0f') '.png']))
end
for kk = 1:size(sETDRS.regionsETDRS_3D{3},3)
    imID = 'ETDRS_3D_3';
    im = sETDRS.regionsETDRS_3D{3}(:,:,kk);
    if ~exist(fullfile(exportFolder,'ExportFiles', folderName, imID), 'dir')
        mkdir(fullfile(exportFolder,'ExportFiles', folderName, imID));
    end
    imwrite(im, fullfile(exportFolder,'ExportFiles', folderName, imID, [imID '_' num2str(kk, '%03.0f') '.png']))
end
for g = 1:size(sETDRS.regionsETDRS_3D{4},3)
    imID = 'ETDRS_3D_4';
    im = sETDRS.regionsETDRS_3D{4}(:,:,g);
    if ~exist(fullfile(exportFolder,'ExportFiles', folderName, imID), 'dir')
        mkdir(fullfile(exportFolder,'ExportFiles', folderName, imID));
    end
    imwrite(im, fullfile(exportFolder,'ExportFiles', folderName, imID, [imID '_' num2str(g, '%03.0f') '.png']))
end
for gg = 1:size(sETDRS.regionsETDRS_3D{5},3)
    imID = 'ETDRS_3D_5';
    im = sETDRS.regionsETDRS_3D{5}(:,:,gg);
    if ~exist(fullfile(exportFolder,'ExportFiles', folderName, imID), 'dir')
        mkdir(fullfile(exportFolder,'ExportFiles', folderName, imID));
    end
    imwrite(im, fullfile(exportFolder,'ExportFiles', folderName, imID, [imID '_' num2str(gg, '%03.0f') '.png']))
end

end

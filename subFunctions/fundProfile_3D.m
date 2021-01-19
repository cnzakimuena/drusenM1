
function volumeProfile = fundProfile_3D(volBW, sETDRS, sizeRed)

% RESULTS
%total drusen volume
% 'radiusFac' is the conversion factor, 3000/1536 or 1.95 um/px 
radiusFac = 3000/1536*1/sizeRed; % um/px
totalVolume = sum(volBW(:))*radiusFac^3/1000^3;
vGR = zeros(1, size(sETDRS.regionsETDRS_3D, 2));
for k = 1:size(sETDRS.regionsETDRS_3D, 2)
    currVol = volBW(sETDRS.regionsETDRS_3D{k});
    realVol = sum(currVol(:))*radiusFac^3/1000^3;
    vGR(:, k) = realVol; % at given grid region, mm^3
end
volumeProfile = [totalVolume vGR]; % mm^3
% figure;imshow3D(volBW,[])

end


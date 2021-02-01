function [DAPIseg] = NuclearWatershed(f,DAPI,Res,DAPIseg)

I_DAPI_bg = imopen(DAPI,strel('disk',(round(Res/10))));
I_DAPI2 = DAPI-I_DAPI_bg;
I_DAPI2_BW = imfill(imbinarize(imadjust(I_DAPI2)),'holes');
I_DAPI2_BW = bwareaopen(I_DAPI2_BW,20);

DAPI_DistTransf1 = bwdist(~I_DAPI2_BW);
DAPI_DistTransf1 = imcomplement(DAPI_DistTransf1);
DAPI_DistTransfMask = imextendedmin(DAPI_DistTransf1,2);
DAPI_DistTransf2 = imimposemin(DAPI_DistTransf1,DAPI_DistTransfMask);
DAPI_Watershed = watershed(DAPI_DistTransf2);
DAPI_Watershed_BW = I_DAPI2_BW;
DAPI_Watershed_BW(DAPI_Watershed == 0) = 0;
DAPI_Watershed_BW2 = logical(imdilate(imerode(DAPI_Watershed_BW,strel('disk',1)),strel('disk',1)));
DAPI_Watershed_BW2 = bwareafilt(DAPI_Watershed_BW2,[20 5000]);

DAPIseg(f).watershedperim = bwperim(imdilate(DAPI_Watershed_BW2,strel('disk',1)));
DAPIseg(f).cc = bwconncomp(DAPI_Watershed_BW2,4);
DAPIseg(f).props = regionprops(DAPIseg(f).cc,DAPI,'MeanIntensity','PixelValues','Area','Centroid','Circularity');

for p = 1:length(DAPIseg(f).props)
DAPIseg(f).NuclearAreas(p,1) = DAPIseg(f).props(p).Area;
end

end
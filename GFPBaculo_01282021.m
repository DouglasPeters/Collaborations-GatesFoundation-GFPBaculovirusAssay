clear all; clc; tic; cd(userpath);

Folder = 'F:\DATA\01.27-01.29.2021 IXM For Gates\01.29.21 IXM Imaging of GFP-Baculo (72hr)\6454\';

Channels = 2;
Ch_DAPI = 1;

FigShow = 1;
FigSave = 0;

%% File sorting and information collection %%
cd(Folder);
srcFilesTotal = dir('*.tif');

START = 1;
FINISH = length(srcFilesTotal)/Channels;

disp('Collecting file information and sorting image channels...');
for i = START:length(srcFilesTotal)
if numel(strfind(srcFilesTotal(i).name,'w1')) > 0
    SortIdx.FindCh1(i,1) = 1;
else SortIdx.FindCh1(i,1) = 0;
end
if numel(strfind(srcFilesTotal(i).name,'w2')) > 0
    SortIdx.FindCh2(i,1) = 1;
else SortIdx.FindCh2(i,1) = 0;
end
if Channels > 2
    if numel(strfind(srcFilesTotal(i).name,'w3')) > 0
        SortIdx.FindCh3(i,1) = 1;
    else SortIdx.FindCh3(i,1) = 0;
    end
else end
if Channels > 3
    if numel(strfind(srcFilesTotal(i).name,'w4')) > 0
        SortIdx.FindCh4(i,1) = 1;
    else SortIdx.FindCh4(i,1) = 0;
    end
else end
end
srcFilesSorted.Ch1 = srcFilesTotal(logical(SortIdx.FindCh1));
srcFilesSorted.Ch2 = srcFilesTotal(logical(SortIdx.FindCh2));
if Channels > 2, srcFilesSorted.Ch3 = srcFilesTotal(logical(SortIdx.FindCh3)); else end
if Channels > 3, srcFilesSorted.Ch4 = srcFilesTotal(logical(SortIdx.FindCh4)); else end
for f = 1:size(srcFilesSorted.Ch1,1) 
    FileInfo.Ch1(f,1).Filename = srcFilesSorted.Ch1(f,1).name; 
    namesort(f,1).underscoreIdx = strfind(FileInfo.Ch1(f,1).Filename,'_');
    namesort(f,1).well = namesort(f,1).underscoreIdx(1,1)+1;
    namesort(f,1).site = namesort(f,1).underscoreIdx(1,2)+2;
    namesort(f,1).channel = namesort(f,1).underscoreIdx(1,3)+2;
    FileInfo.Ch1(f,1).Well = FileInfo.Ch1(f,1).Filename(namesort(f,1).well:(namesort(f,1).well+2));
    FileInfo.Ch1(f,1).Site = FileInfo.Ch1(f,1).Filename(namesort(f,1).site:namesort(f,1).site);
    FileInfo.Ch1(f,1).Channel = FileInfo.Ch1(f,1).Filename(namesort(f,1).channel:namesort(f,1).channel); 
    FileInfo.Ch2(f,1).Filename = srcFilesSorted.Ch2(f,1).name;
    namesort(f,1).underscoreIdx = strfind(FileInfo.Ch2(f,1).Filename,'_');
    namesort(f,1).well = namesort(f,1).underscoreIdx(1,1)+1;
    namesort(f,1).site = namesort(f,1).underscoreIdx(1,2)+2;
    namesort(f,1).channel = namesort(f,1).underscoreIdx(1,3)+2;
    FileInfo.Ch2(f,1).Well = FileInfo.Ch2(f,1).Filename(namesort(f,1).well:(namesort(f,1).well+2));
    FileInfo.Ch2(f,1).Site = FileInfo.Ch2(f,1).Filename(namesort(f,1).site:namesort(f,1).site);
    FileInfo.Ch2(f,1).Channel = FileInfo.Ch2(f,1).Filename(namesort(f,1).channel:namesort(f,1).channel);  
    
    SiteNum(f,1) = str2num(FileInfo.Ch1(f,1).Site); 
end
    SiteMax = max(SiteNum);

if FigShow == 1, figure; else end
for f = START:FINISH
        
    time(f,1).ElapsedSeconds = toc;
    
   try
    clc
    progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2);
    
    if f == START
        cd(Folder); mkdir('Analysis'); cd(Folder);
    else end

    if progress < 10
        disp('Estimated time remaining will display after 10% of images are analyzed...');
    else
        time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/((FINISH-START+1)-(FINISH-f));
        time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH-START+1);
        time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
        time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
        time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
        estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
        estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
        disp(estimate);
        disp(estimate2);
    end
   
    
    I.Ch1 = imread(FileInfo.Ch1(f,1).Filename);
    I.Ch2 = imread(FileInfo.Ch2(f,1).Filename);
    if Channels > 2, I.Ch3 = imread(FileInfo.Ch3(f,1).Filename); else end
    if Channels > 3, I.Ch4 = imread(FileInfo.Ch4(f,1).Filename); else end
    
    Res = size(I.Ch1,1);
    
    %% Nuclear Segmentation and Watershedding %%
    disp('Segmenting and Watershedding Nuclei...');
    if Ch_DAPI == 1, DAPI = I.Ch1;
    elseif Ch_DAPI == 2, DAPI = I.Ch2;
    elseif Ch_DAPI == 3, DAPI = I.Ch3;
    elseif Ch_DAPI == 4, DAPI = I.Ch4;
    else warning('Invalid "Ch_DAPI" value.');
    end
    
    if f==START, DAPIseg = struct(); else end    
    [DAPIseg] = NuclearWatershed(f,DAPI,Res,DAPIseg);
    
    %% Cellular Analysis %%
    disp('Extracting cellular properties...');
    
    Results(f,1).Filename = FileInfo.Ch1(f).Filename;
    Results(f,1).Well = FileInfo.Ch1(f).Well;
    Results(f,1).Site = FileInfo.Ch1(f).Site;
    Results(f,1).NumberNuclei = DAPIseg(f).cc.NumObjects;
    
    Results(f,1).Ch1Props = regionprops(DAPIseg(f).cc,I.Ch1,'Area','MeanIntensity');
    Results(f,1).Ch2Props = regionprops(DAPIseg(f).cc,I.Ch2,'Area','MeanIntensity');
    if Channels > 2, Results(f,1).Ch3Props = regionprops(DAPIseg(f).cc,I.Ch3,'Area','MeanIntensity'); else end
    if Channels > 3, Results(f,1).Ch4Props = regionprops(DAPIseg(f).cc,I.Ch4,'Area','MeanIntensity'); else end
    
    for m = 1:size(Results(f).NumberNuclei)
    Results(f,1).Ch1MeanInt(m,1) = Results(f).Ch1Props(m).MeanIntensity;
    Results(f,1).Ch2MeanInt(m,1) = Results(f).Ch2Props(m).MeanIntensity;
    if Channels > 2, Results(f,1).Ch3MeanInt(m,1) = Results(f).Ch3Props(m).MeanIntensity; else end
    if Channels > 3, Results(f,1).Ch4MeanInt(m,1) = Results(f).Ch4Props(m).MeanIntensity; else end
    end
        

    %% Figure %%
    if FigShow == 1
    disp('Generating Figure...');
    
    figimage1 = imoverlay(imadjust(DAPI),imdilate(DAPIseg(f).watershedperim,strel('disk',1)),[0.7 0.7 0]);
    figimage2 = imoverlay(imadjust(I.Ch2),imdilate(DAPIseg(f).watershedperim,strel('disk',1)),[0.7 0.7 0]);
    
    C = [figimage1 figimage2];
    imshow(C); title('Segmentation Summary');
    else end


    catch
end
end 

%% Results Sorting %%
% for r = 1:numel(SiteNum)
%     for s = 1:SiteMax
%         SiteIdx(r,s) = SiteNum(r) == s;
%     end
% end

disp('Saving Results Files...');
cd(Folder); cd Analysis;
save('AnalysisResults.mat','Results','-v7.3');



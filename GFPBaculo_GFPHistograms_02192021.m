clear all; clc; tic; cd(userpath);

Folder = 'F:\DATA\Gates02192021\2021-02-19\6484\';

Ch_GFP = 2;
Channels = 2;

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

for f = START:FINISH
        
    time(f,1).ElapsedSeconds = toc;
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
    
    if Ch_GFP == 1, GFP_FileInfo = FileInfo.Ch1(f,1); else end
    if Ch_GFP == 2, GFP_FileInfo = FileInfo.Ch2(f,1); else end
    if Ch_GFP == 3, GFP_FileInfo = FileInfo.Ch3(f,1); else end
    if Ch_GFP == 4, GFP_FileInfo = FileInfo.Ch4(f,1); else end
    
    GFP = imread(GFP_FileInfo.Filename);
    Res = size(GFP,1);
    
    
    Results(f,1).Filename = GFP_FileInfo.Filename;
    Results(f,1).Well = GFP_FileInfo.Well;
    Results(f,1).Site = str2num(GFP_FileInfo.Site);
  
    %% GFP Signal Distribution Analysis %%

    if Results(f,1).Site == 1, GFP_Merge(:,:,1) = GFP;
    elseif Results(f,1).Site == 2, GFP_Merge(:,:,2) = GFP;
    elseif Results(f,1).Site == 3, GFP_Merge(:,:,3) = GFP;
    elseif Results(f,1).Site == 4, GFP_Merge(:,:,4) = GFP;
    elseif Results(f,1).Site == 5
        GFP_Merge(:,:,5) = GFP;
        GFP_Histo = histogram(GFP,'Normalization','Probability');
        GFP_Histo.BinWidth = 5000;
        GFP_Histo.BinLimits = [5000,40000];
        Results(f,1).GFPDistribution(1:length(GFP_Histo.Values),1) = GFP_Histo.Values;
        drawnow;
    else warning('Assuming there are only 5 sites per well....');
    end

end 
for w = 1:(length(Results)/SiteMax)
    Sorted(:,w) = Results(w*5,1).GFPDistribution;
end


disp('Saving Results Files...');
cd(Folder); cd Analysis;
save('AnalysisResults.mat','Results','-v7.3');



%% Single Cell Analysis
% Conversion factor: 73 fluo light intensity = 1 mRNA

clear all
close all
clc


numFrames = 121;    % 4 hours
samplingTime = 2;   % minutes
reference = 1500;
timeAxis = [0:(numFrames-1)]*samplingTime;
refVec = ones(1,numFrames)*reference;
colorsExp = [55,126,184
    228,26,28
    245,170,66
    221,245,66
    141,245,66
    66,245,170
    ]/255;

expNum = [45];

for loopNum=1:length(expNum)
    
    load(['data' num2str(expNum(loopNum)) '.mat']);
    
    %% Extact only the relevant cells (which are segmented in all frames)
    relevantCellIndx = [];
    relevantIndxCount = 0;
    for i=1:length(cellTrajectory)
        if length(cellTrajectory{i}.spotIntensity)==numFrames
            checkMinusOne = find(cellTrajectory{i}.spotIntensity==-1);
            if isempty(checkMinusOne)
                relevantIndxCount = relevantIndxCount + 1;
                relevantCellIndx(relevantIndxCount) = i;
            end
        end
    end
    relevantCellTrajectory = cellTrajectory(relevantCellIndx);
    
    %% Cell Intensity and Control Output
    cellIntensity = zeros(relevantIndxCount, numFrames);
    controlOutput = zeros(relevantIndxCount, numFrames);
    inputIntensity = zeros(relevantIndxCount, numFrames);
    err = zeros(relevantIndxCount, numFrames);
    controlV = zeros(relevantIndxCount, numFrames);
    for i=1:relevantIndxCount
        cellIntensity(i,:) = relevantCellTrajectory{i}.spotIntensity;
        inputIntensity(i,:) = relevantCellTrajectory{i}.inputIntensity;
    end
    meanCellIntensity = mean(cellIntensity);
    meanInputIntensity = mean(inputIntensity);
    
    %%
    fig1 = figure(1);
    hold on;
    
    fltrWindow = 5;
    meanCellIntensityAvg = conv(meanCellIntensity,ones(1,fltrWindow)/fltrWindow,'valid');
    p(loopNum) = plot(timeAxis((1+(fltrWindow-1)/2):(end-(fltrWindow-1)/2)), meanCellIntensityAvg/73, 'Color', colorsExp(loopNum,:),'LineWidth',2);
    plot(timeAxis, meanCellIntensity/73, 'Color', colorsExp(loopNum,:), 'LineStyle', ':','LineWidth',1);
    title('Mean output')
    xlabel('time [minutes]')
    ylabel('# nascent mRNAs')
    ylim([0 3400/73])
    xlim([0 240])
    hold off;
        
end

figure(1); hold on;
ref = refline(0,reference/73);
ref.LineStyle = '--';
ref.LineWidth = 2;
ref.Color = 'k';
hold off;

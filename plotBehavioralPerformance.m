function plotBehavioralPerformance()


person = 'EEL';
mouse = [person, '13'];
pth = ['Z:\' person '\Experiments\' mouse '\Analysis'];
dates = {'2020-09-09'};
ix = 1;

params.dates = dates; 
params.datesIX = 1;
% params.tag = 'Spontaneous Water';
params.tag = '2AFC Autowater';


% obj = getDataStructures(pth);
obj = load(fullfile(pth,dates{1},'bpod.mat'));
% dates = {obj.date};



% biParams.dates = dates(1:3)'; % specify desired dates, can input indices of dates or dates in form YEAR-MO-DA
% biParams.datesIX = getDateIX(dates, biParams); % get index of desired dates
% biParams.tag = 'Bilateral';
% 
% uniParams.dates = dates(4:5)';
% uniParams.datesIX = getDateIX(dates, uniParams);
% uniParams.tag = 'Unilateral';



redColor = [1 0.05 0.05];
blueColor = [0 0.5 0.95];

% whichTT = {
%     'tt.hit', 'tt.L&~tt.no&tt.stim.num==0', 1, redColor;
%     'tt.hit', 'tt.R&~tt.no&tt.stim.num==0', 1, blueColor;
%     'tt.hit', 'tt.L&~tt.no&tt.stim.num==1', 2, redColor;
%     'tt.hit', 'tt.R&~tt.no&tt.stim.num==1', 2, blueColor;
%     'tt.hit', 'tt.L&~tt.no&tt.stim.num==2', 3, redColor;
%     'tt.hit', 'tt.R&~tt.no&tt.stim.num==2', 3, blueColor;
% %     'tt.hit', 'tt.L&~tt.no&tt.stim.num==3', 4, redColor;
% %     'tt.hit', 'tt.R&~tt.no&tt.stim.num==3', 4, blueColor;
% %     'tt.hit', 'tt.L&~tt.no&tt.stim.num==4', 5, redColor;
% %     'tt.hit', 'tt.R&~tt.no&tt.stim.num==4', 5, blueColor
%     };

% biParams.redColor = redColor;
% biParams.blueColor = blueColor;
% biParams.whichTT = whichTT;
% uniParams.redColor = redColor;
% uniParams.blueColor = blueColor;
% uniParams.whichTT = whichTT;

% perf = getBehavioralPerf(obj, params);

% biPerf  = getBehavioralPerf(obj, biParams);
% uniPerf  = getBehavioralPerf(obj, uniParams);


% bBil = makeBarGraph(biPerf, biParams, mouse);
% bUni = makeBarGraph(uniPerf, uniParams, mouse);

% plot lick rasters for all days
for i=1:numel(obj)
    ctrlTs = find(obj.bp.autowater.nums==2);
    autoTs = find(obj.bp.autowater.nums==1);
%     blocks = {1:399, 400:obj(i).bp.Ntrials};
    blocks = {ctrlTs, autoTs};
    plotLickRasters(obj(i).bp, redColor, blueColor, mouse, dates{i}, params.tag, blocks)
    plotLick1Hists(obj(i).bp, redColor, blueColor, mouse, dates{i}, params.tag, blocks)
end



% plotBlockLickRasters(obj(i).bp, redColor, blueColor, mouse, dates{i}, params.tag)

% plot rasters for all days in params.dates
% plotSummaryRasters(biParams, obj, mouse)
% plotSummaryRasters(uniParams, obj, mouse)


function plotLick1Hists(bp, redColor, blueColor, mouse, date, tag, blocks)
params.redColor = redColor;
params.blueColor = blueColor;


ctrlTrials = blocks{1};
if length(blocks)==2
    expTrials = blocks{2};
end

numStims = numel(unique(bp.stim.num));

figure;
for i=1:numStims
    subplot(numStims,2,2*i-1); hold on;
    figTitle = ['Left Trials Stim ', num2str(i-1)];
    plotLickHists(bp, find(bp.L&bp.stim.num==(i-1)), params, figTitle, ctrlTrials);
    subplot(numStims,2,2*i); hold on;
    figTitle = ['Right Trials Stim ', num2str(i-1)];
    plotLickHists(bp, find(bp.R&bp.stim.num==(i-1)), params, figTitle, ctrlTrials);
     
    sgtitle(['First Licks, ', mouse, ' ', date, ' 2AFC'], 'FontWeight', 'bold', 'FontSize', 14);
end

if exist('expTrials', 'var')
    figure;
    for i=1:numStims
        subplot(numStims,2,2*i-1); hold on;
        figTitle = ['Left Trials Stim ', num2str(i-1)];
        plotLickHists(bp, find(bp.L&bp.stim.num==(i-1)), params, figTitle, expTrials);
        subplot(numStims,2,2*i); hold on;
        figTitle = ['Right Trials Stim ', num2str(i-1)];
        plotLickHists(bp, find(bp.R&bp.stim.num==(i-1)), params, figTitle, expTrials);
        
        sgtitle(['First Licks, ', mouse, ' ', date, ' ', tag], 'FontWeight', 'bold', 'FontSize', 14);
    end
end


function plotLickHists(bp, trials, params, figureTitle, blocks)
trials = trials(ismember(trials,blocks));

lick1L = NaN*zeros(1,length(trials));
lick1R = NaN*zeros(1,length(trials));
for i=1:length(trials)
    
    if ~isempty(bp.ev.lickL{trials(i)})
        lick1L(i) = bp.ev.lickL{trials(i)}(1)-bp.ev.goCue(trials(i));
    end
    if ~isempty(bp.ev.lickR{trials(i)})
        lick1R(i) = bp.ev.lickR{trials(i)}(1)-bp.ev.goCue(trials(i));
    end
    
    if ~isnan(lick1L(i)) && ~isnan(lick1R(i)) %if licked to both sides on trial
        if lick1L(i) < lick1R(i) % only keep first lick
            lick1R(i) = NaN;
        else
            lick1L(i) = NaN;
        end
    end
    
end
histogram(lick1L, 'BinWidth', 0.005, 'FaceColor', params.redColor); hold on;
histogram(lick1R, 'BinWidth', 0.005, 'FaceColor', params.blueColor);
title(figureTitle, 'FontSize', 12, 'FontWeight', 'Light')
xlim([0.05 0.25])



function plotSummaryRasters(params, obj, mouse)
dates = params.dates;
datesIX = params.datesIX;

figh = figure;
set(figh, 'Position', [1000 200 200*numel(dates) 800])
for i=1:numel(dates)
    bp = obj(datesIX(i)).bp;
    
    subplot(2, numel(dates), i); hold on;
    plotByLick(bp, find(bp.L), params, '', dates(i));
    title(['L ' obj(datesIX(i)).date(end-4:end)], 'FontWeight', 'normal')
    set(gca, 'FontSize', 10)
    
    subplot(2, numel(dates), i+numel(dates)); hold on;
    plotByLick(bp, find(bp.R), params, '', dates(i));
    title(['R ' obj(datesIX(i)).date(end-4:end)], 'FontWeight', 'normal')
    set(gca, 'FontSize', 10)
end
sgtitle([mouse, ' ', params.tag], 'FontWeight', 'bold')





function plotLickRasters(bp, redColor, blueColor, mouse, date, tag, blocks)
params.redColor = redColor;
params.blueColor = blueColor;


ctrlTrials = blocks{1};
if length(blocks)==2
    expTrials = blocks{2};
end

figure;
subplot(1,2,1); hold on;
plotByLick(bp, find(bp.L), params, 'Left Trials', date, ctrlTrials);

subplot(1,2,2); hold on;
plotByLick(bp, find(bp.R), params, 'Right Trials', date, ctrlTrials);

sgtitle([mouse, ' ', date, ' 2AFC'], 'FontWeight', 'bold'); set(gca, 'FontSize', 14)

if exist('expTrials', 'var')
    if ~isempty(expTrials)
        figure;
        subplot(1,2,1); hold on;
        plotByLick(bp, find(bp.L), params, 'Left Trials', date, expTrials);
        
        subplot(1,2,2); hold on;
        plotByLick(bp, find(bp.R), params, 'Right Trials', date, expTrials);
        
        sgtitle([mouse, ' ', date, ' ', tag], 'FontWeight', 'bold'); set(gca, 'FontSize', 14)
    end
end





function plotByLick(bp, alltrials, params, figureTitle, date, blocks)
% figure; hold on;
stimtype = unique(bp.stim.num);
Nplotted = 0;
for i = 1:numel(stimtype)
    
    trials = alltrials(bp.stim.num(alltrials)==stimtype(i));
    if i==1
        trials = alltrials(bp.stim.num(alltrials)==stimtype(i));
    elseif i==2
        trials = alltrials(bp.stim.num(alltrials)==stimtype(i)&strcmp(bp.stim.state(alltrials),'Delay'));
    end
    trials = trials(ismember(trials,blocks));
    
    y = [Nplotted Nplotted+numel(trials)];
    
%     if strcmp(date, '2020-02-27')
%         if stimtype(i) == 1
%             stimtype(i) = 3;
%         end
%     end
    
    uniColor = [0.85 0.85 1];
    biColor = [0.7 0.7 1];
    
    if stimtype(i)==1
        myfill([-0.9 -0.2], y+0.5, 0.2, uniColor);
    elseif stimtype(i)==2
        myfill([-0.9 -0.2], y+0.5, 0.2, biColor);
%     elseif stimtype(i)==3
%         myfill([2.5 3.5], y, 0.2);
%     elseif stimtype(i)==4
%         myfill([3 4], y, 0.2);
    end
    
    for j = 1:numel(trials)
        
        lickL = bp.ev.lickL{trials(j)}-bp.ev.goCue(trials(j));
        lickR = bp.ev.lickR{trials(j)}-bp.ev.goCue(trials(j));
        
        plot(lickL, Nplotted+j*ones(size(lickL)), '.', 'Color', params.redColor); hold on;
        plot(lickR, Nplotted+j*ones(size(lickR)), '.', 'Color', params.blueColor);
        
        if bp.early(trials(j))
            plot([0 10], Nplotted+j+[0 0], 'k:');
        end
        
    end
    Nplotted = Nplotted+numel(trials);
%     plot(bp.ev.sample(trials), 1:numel(trials), 'k:');
%     plot(bp.ev.delay(trials), 1:numel(trials), 'k:');
%     plot(bp.ev.goCue(trials), 1:numel(trials), 'k:');

    plot([0 10], Nplotted+[0.5 0.5], 'k-');
    
    if i==2
        trials = alltrials(bp.stim.num(alltrials)==stimtype(i)&strcmp(bp.stim.state(alltrials),'GoCue'));
        trials = trials(ismember(trials,blocks));
        
        y = [Nplotted Nplotted+numel(trials)];
        
        %     if strcmp(date, '2020-02-27')
        %         if stimtype(i) == 1
        %             stimtype(i) = 3;
        %         end
        %     end
        
        uniColor = [0.85 0.85 1];
        biColor = [0.7 0.7 1];
        
        if stimtype(i)==1
            myfill([0 0.7], y+0.5, 0.2, uniColor);
        elseif stimtype(i)==2
            myfill([0 0.7], y+0.5, 0.2, biColor);
            %     elseif stimtype(i)==3
            %         myfill([2.5 3.5], y, 0.2);
            %     elseif stimtype(i)==4
            %         myfill([3 4], y, 0.2);
        end
        
        for j = 1:numel(trials)
            
            lickL = bp.ev.lickL{trials(j)}-bp.ev.goCue(trials(j));
            lickR = bp.ev.lickR{trials(j)}-bp.ev.goCue(trials(j));
            
            plot(lickL, Nplotted+j*ones(size(lickL)), '.', 'Color', params.redColor); hold on;
            plot(lickR, Nplotted+j*ones(size(lickR)), '.', 'Color', params.blueColor);
            
            if bp.early(trials(j))
                plot([0 10], Nplotted+j+[0 0], 'k:');
            end
            
        end
        Nplotted = Nplotted+numel(trials);
        %     plot(bp.ev.sample(trials), 1:numel(trials), 'k:');
        %     plot(bp.ev.delay(trials), 1:numel(trials), 'k:');
        %     plot(bp.ev.goCue(trials), 1:numel(trials), 'k:');
        
        plot([0 10], Nplotted+[0.5 0.5], 'k-');
    end
end

plot([0 10], Nplotted+[0.5 0.5], 'k-', 'LineWidth', 2);

xlabel('Time (s)');
title(figureTitle)
set(gca, 'FontSize', 12);
xlim([-1 6]);
% ylim([0 Nplotted]);


function myfill(xl, yl, rampdur, color)
x = [xl(1) xl(1) xl(2) xl(2)];
y = [yl(1) yl(2) yl(2) yl(1)];
f = fill(x, y, color); set(f, 'LineStyle', 'None');
f = fill([x(3:4) x(3:4)+rampdur], y, [0.85 0.85 1]); set(f, 'LineStyle', 'None');






function datesIX = getDateIX(dates, params)
datesIX = [];
for i=1:numel(dates)
    if ismember(dates{i}, params.dates)
        datesIX = [datesIX, i];
    end
end



function perf = getBehavioralPerf(obj, params)
% whichTT in form {'tt.numerator', 'tt.denominator'; ...}
% if denominator is all trials, leave empty
numTT = size(params.whichTT,1);
numDates = numel(params.dates);
whichTT = params.whichTT;

numTrials = NaN.*ones(numTT,numDates);
denTrials = NaN.*ones(numTT,numDates);
perc = NaN.*ones(numTT,numDates);
sem = NaN.*ones(numTT,numDates);

for i=1:numDates
    tt = obj(params.datesIX(i)).bp;
    for j=1:numTT
        
        numTrials(j,i) = eval(['nansum(' whichTT{j,1} '&' whichTT{j,2} ')']);
        if ~isempty(whichTT(j,2))
            denTrials(j,i) = eval(['nansum(' whichTT{j,2} ')']);
        else
            denTrials(j,i) = eval(['numel(' whichTT{j,1} ')']);
        end
        
        if denTrials(j,i)~=0
            perc(j,i) = numTrials(j,i)/denTrials(j,i);
            ssd = sqrt(((1-perc(j,i))^2.*numTrials(j,i) + perc(j,i)^2.*(denTrials(j,i)-numTrials(j,i)))./denTrials(j,i));
            sem(j,i) = ssd/sqrt(denTrials(j,i));
        else
            perc(j,i) = NaN;
        end
    end
end

perf.numTrials = numTrials;
perf.denTrials = denTrials;
perf.perc = perc;
perf.sem = sem;




function barObj = makeBarGraph(perf, params, mouse)
figure;
% bar graph with individual date performances and standard error of the
% mean overlaid
whichTT = params.whichTT;
perc = perf.perc;
sem = perf.sem;

xpoints = [];
% manual grouping for bar chart
for i=1:size(whichTT,1)
    if mod(i,2)==1
        xpoints = [xpoints i+0.1];
    else
        xpoints = [xpoints i-0.1];
    end
end
barObj = bar(xpoints, nanmean(perc,2));
barObj.FaceColor = 'flat';
barObj.CData(:,:) = reshape([whichTT{:,4}], 3, size(whichTT,1))';
barObj.EdgeColor = 'none';
ax1 = gca;
ylabel('Performance')
title([mouse, ' ', params.tag])
set(ax1, 'ylim', [0 1], 'XTick', [1.5, 3.5, 5.5, 7.5, 9.5], 'XTickLabel', ...
    {'Control', 'Stim 1', 'Stim 2', 'Stim 3', 'Stim 4'}, 'FontSize', 14,...
    'YGrid', 'on', 'YMinorGrid', 'on')
hold(ax1, 'all')

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','none','YColor','none',...
    'XTick', [], 'YTick', [], 'xlim', get(ax1, 'xlim'), 'ylim', get(ax1, 'ylim'));
hold(ax2, 'all')

for i=1:size(whichTT,1)
    numDataPts = sum(~isnan(perc(i,:)));
    percX = xpoints(i).*ones(1,numDataPts)+linspace(-0.1,0.1,numDataPts); %linspace keeps data points from overlapping
    sz = 50;
    scatter(percX,perc(i,~isnan(perc(i,:))), sz, 'Parent', ax2, 'k', 'filled')
%     for j=1:numel(percX)
%         errorbar(percX,perc(i,~isnan(perc(i,:))),sem(i,~isnan(perc(i,:))), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1)
%     end
end

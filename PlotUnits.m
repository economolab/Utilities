function PlotUnits()
set(0, 'DefaultFigurePosition', [1138          84        1392        1245])

pth = 'Y:\EEL\Experiments\EEL14\Analysis\';

obj = getDataStructures(pth);


% plotPSTHs(obj(1),1)


% load(fullfile('E:\EEL6\Analysis\2020-03-01', 'data_structure_EEL6_2020-03-01.mat'));
% % load(fullfile('Y:\EEL\Experiments\EEL2\Analysis\2019-09-12', 'data_structure_EEL2_2019-09-12.mat'));

probe = 1;
% plotPSTHs(obj(2), probe)

for i = 1:numel(obj)
    thisobj = obj(i);
    R = thisobj.bp.R;
    L = thisobj.bp.L;
    hit = thisobj.bp.hit;
    early = thisobj.bp.early;
    lickL = thisobj.bp.ev.lickL;
    lickR = thisobj.bp.ev.lickR;
    stim = thisobj.bp.stim;
    no = thisobj.bp.no;
    
    med = 7;
    view =2;
    bodyparts = 1:3;
    offset = 25;
    figure;
    plotStackedTrajectories(thisobj, bodyparts(1), view, {L&stim.num==1; L&stim.num==2; L&stim.num==3; L&stim.num==4}, med, {[1 0 0]; [0 0 1]; [0 0 0]; [0 0.75 0]}, offset);
end



figure;
plotTrajectories(thisobj, bodyparts, view, R&hit, med, 'b');
plotTrajectories(thisobj, bodyparts, view, R&no, med, 'r');







function plotProjections(obj)
probe = 1;
dt = 0.005;
edges = -3:dt:3;
Nsess = numel(obj);

tt = cell(Nsess, 1);
psth = cell(Nsess, 1);
meta = cell(Nsess, 1);
ob = cell(Nsess, 1);

for i = 3
    

    [tt{i}, psth{i}, meta{i}, ob{i}] = ObjToPSTH(obj(i), probe, edges);
%     tt{i}.filename = fn{i};
%     meta{i}.filename = fn{i};
%     meta{i}.animal = animal;

end

psth = psth{1};
tt = tt{1};

goodclu = mean(mean(psth, 3), 1)>0.1;
psth = psth(:, goodclu, :);








function [tt, psth, meta, obj] = ObjToPSTH(obj, probe, edges)

dt = mean(diff(edges));

Nunits = numel(obj.clu{probe});


quality = cell(Nunits, 1);
channel = zeros(Nunits, 1);
for j = 1:Nunits
%     quality{j} = unit.quality;
    if isempty(obj.clu{probe}(j).site)
        channel(j) = 1;
    else
        channel(j) = obj.clu{probe}(1).site;
    end
end

tt = obj.bp;
tt.lickL = obj.bp.R&obj.bp.miss | obj.bp.L&obj.bp.hit;
tt.lickR = obj.bp.L&obj.bp.miss | obj.bp.R&obj.bp.hit;


tt.trialtypenames = unique(obj.bp.stim.num);
NtrialTypes = numel(tt.trialtypenames);
tt.trialtypeflags = cell(NtrialTypes, 1);

for i = 1:NtrialTypes

    tt.trialtypeflags{i} = obj.bp.stim.num==i;
    
end


meta.trialTypeStr = tt.trialtypenames;
% meta.trialOutcomeStr = obj.trialTypeStr;
% meta.trialOutcomeFlags = obj.trialTypeMat;
% meta.unitQuality = quality;
meta.channel = channel;
% siteMap = getSiteLocs(obj.sessionMeta.probeType, obj.sessionMeta.depth);
% meta.depth = siteMap(meta.channel);
% meta.manipulatorDepth = obj.sessionMeta.depth;
meta.unitNum = numel(obj.clu{probe});
Ntrials = obj.bp.Ntrials;

psth = NaN.*zeros(numel(edges), Nunits, Ntrials);

cnt = 0;
meta.unitNumber = Nunits;
for j = 1:Nunits
    cnt = cnt+1;
    
    unit = obj.clu{probe}(j);
    
    spktrial = unit.trial;
    spktrialtm = unit.trialtm;
    
    %     trial1 = max(min(spktrial), find(tt.hit, 1, 'first')+10);
    %     trial2 = min(max(spktrial), find(tt.hit, 1, 'last')-10);
    
    trial1 = 1;
    trial2 = Ntrials;
    
    for i = trial1:trial2
        ix = spktrial==i;%obj.trialIDs(i);
        cuetm = obj.sglx.rewardFileOffset;
        t = spktrialtm(ix)-cuetm(i);
        
        psth(:, cnt, i) = histc(t, edges)./dt;
        
    end
    
end
psth = psth(1:end-1, :, :);









 
function plotStackedTrajectories(obj, bodypart, view, trials, med, clrs, offset)

Ntrials = numel(obj.traj{1});

dt = 1./400;
yoffset = 0;
for k = 1:2
    subplot(1,2,k); hold on;
    for j = 1:numel(trials)
        for i = 1:Ntrials
            if ~any(trials{j}(i))
                continue;
            end
            
            ts = medfilt1(obj.traj{view}(i).ts(:, k, bodypart), med, [], 1);
            time = dt.*(1:numel(ts));
            plot(time, ts+yoffset*offset, 'Color', clrs{j});
            %     xoffset = xoffset+time(end);
            yoffset = yoffset+1;
        end
    end
end
% xlim([0 5]);
yl = ylim();
plot(2.48+[0 0], yl, 'k:');



% function out = getDataStructures(pth)
% data_dirs = dir(pth);
% for i=1:numel(data_dirs)
%    if length(data_dirs(i).name) == 10
%        data_files = dir(fullfile(pth, data_dirs(i).name));
%        for j = 1:numel(data_files)
%            if contains(data_files(j).name, 'data_structure')
%                
%               fn{i} = fullfile(pth, data_dirs(i).name, data_files(j).name);
% 
%            end
%        end
%    end
% end
% fn = fn(~cellfun(@isempty, fn));
% 
% for i=1:numel(fn)
%    a(i) = load(fn{i}); 
% end
% 
% for i = 1:numel(a)
%     out(i) = a(i).obj;
% end




function plotPSTHs(obj,probe)

Nclu = numel(obj.clu{probe});
R = obj.bp.R;
L = obj.bp.L;
hit = obj.bp.hit;
early = obj.bp.early;
lickL = obj.bp.ev.lickL;
lickR = obj.bp.ev.lickR;
stim = obj.bp.stim;


for i = 1:Nclu
    c = obj.clu{probe}(i);
    
    [Rix, Rnums] = ismember(c.trial, find(R&hit));
    [Lix, Lnums] = ismember(c.trial, find(L&hit));
    
    [REix, REnums] = ismember(c.trial, find(R&~hit));
    [LEix, LEnums] = ismember(c.trial, find(L&~hit));
    
    binedges = -10:0.05:10;
    binmid = binedges+mean(diff(binedges))/2;
    
    NR = histc(c.trialtm(Rix), binedges);
    NL = histc(c.trialtm(Lix), binedges);
    
    figure(1);
    subplot(3,1,1);
    plot(c.trialtm(Rix), Rnums(Rix), 'b.'); xlim([0 6]);
    subplot(3,1,2);
    plot(c.trialtm(Lix), Lnums(Lix), 'r.');xlim([0 6]);
    subplot(3,1,3);
    plot(binmid, NR, 'b', binmid, NL, 'r'); xlim([0 6]);
pause(0.2);
end

% figure; plot(obj.clu{probe}(i).trialtm, obj.clu{probe}(i).trial, 'r.')
% hold on;
% for k = 1:max(obj.clu{probe}(i).trial)
%     
%     lickR = obj.bp.ev.lickR{k};
%     lickL = obj.bp.ev.lickL{k};
% 
%     
%     plot(lickR, k*ones(size(lickR)), 'ko');
%     plot(lickL, k*ones(size(lickL)), 'ko');
%     
% end





function ts = getTrajectories(obj, bodypart, view, trials)

Ntrials = sum(trials);
len = 0;
for i = 1:numel(obj.traj{view})
    len = max(len, size(obj.traj{view}(i).ts, 1));
end

ts = NaN.*zeros(len, 3, Ntrials);
cnt = 0;

trialnum = find(trials);
for i = 1:Ntrials
    
    cnt = cnt+1;
    dat = obj.traj{view}(trialnum(i)).ts(:, :, bodypart);
    ts(1:size(dat, 1), :, cnt) = dat;
end



function ts = plotTrajectories(obj, bodyparts, view, trials, med, clr)


for i = bodyparts
    bodypart = i;
    
    hold on;
    
    ts = medfilt1(getTrajectories(obj, bodypart, view, trials), med, [], 1);
    plot(squeeze(ts(:, 1, :)), -squeeze(ts(:, 2, :)), 'Color', clr)

end



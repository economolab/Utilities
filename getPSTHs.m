function [tt, psth, meta, ob] = getPSTHs()

dt = 0.005;
parent ='E:\EEL6\Analysis\2020-02-27';
fn = {'data_structure_EEL6_2020-02-27.mat'};
edges = -3:dt:3;

tt = cell(numel(fn), 1);
psth = cell(numel(fn), 1);
meta = cell(numel(fn), 1);
ob = cell(numel(fn), 1);

for i = 1:numel(fn)
    
    load(fullfile(parent, fn{i}));
    [tt{i}, psth{i}, meta{i}, ob{i}] = ObjToPSTH(obj, edges);
    tt{i}.filename = fn{i};
    meta{i}.filename = fn{i};
%     meta{i}.animal = animal;

end

psth = psth{1};
tt = tt{1};
meta = meta{1};


params.ix1 = 300:400;
params.ix2 = 500:550;

params.trial1 = tt.R&tt.hit&tt.stim.num==0;
params.trial2 = tt.L&tt.hit&tt.stim.num==0;

params.probe = [1];

cv = calcCodingVector(psth, meta, params);

trialsToPlot = 1:275;
figure; hold on;
plotCVProj(obj, psth, cv, tt.L&tt.hit&tt.stim.num==0, [1 0 0], trialsToPlot, 50);
plotCVProj(obj, psth, cv, tt.R&tt.hit&tt.stim.num==0, [0 0 1], trialsToPlot, 50 );
% plotCVProj(psth, cv, tt.R&tt.miss&tt.stim.num==0, [0 1 1], trialsToPlot);

plotCVProj(psth, cv, tt.R&tt.miss, [0 1 1], trialsToPlot);
plotCVProj(psth, cv, tt.L&tt.miss, [1 0 1], trialsToPlot);   


plotWeightedNeurons(psth, cv, tt)







function PCAstuff(psth, tt)

L = MySmooth(mean(psth(:,:, tt.L&tt.hit),3),40);
R = MySmooth(mean(psth(:,:, tt.R&tt.hit),3),40);

ix = 10:60;
for i=1:size(psth, 2)
    Lz(:,i) = (L(:,i)-mean(L(ix,i)))/std(L(ix,i));
    Rz(:,i) = (R(:,i)-mean(R(ix,i)))/std(R(ix,i));
end

Lz(isnan(Lz)) = 0;
Rz(isnan(Rz)) = 0;

Lz(isinf(Lz)) = 0;
Rz(isinf(Rz)) = 0;

dat = [Lz; Rz];
dat = dat - repmat(mean(dat, 1), size(dat, 1), 1);

[COEFF, SCORE] = pca(dat, 'NumComponents', 10);

figure; hold on;
plotPCAProj(psth, COEFF, tt.L&tt.hit, 'r');
plotPCAProj(psth, COEFF, tt.R&tt.hit, 'b');

plotPCAProj(psth, COEFF, tt.L&tt.stim.num==3, 'm');
plotPCAProj(psth, COEFF, tt.R&tt.stim.num==3, 'c');



function plotWeightedNeurons(psth, cv, tt)
[~, ix] = sort(abs(cv), 'descend');


figure(1);
for i=1:25
    unit = ix(i);
    ts = squeeze(psth(:,unit,:));
    subplot(5,5,i)
    hold on;
    plot(MySmooth(mean(ts(:,tt.R&tt.hit&tt.stim.num==0),2),100))
    plot(MySmooth(mean(ts(:,tt.R&tt.hit&tt.stim.num==1),2),100))
    plot(MySmooth(mean(ts(:,tt.R&tt.hit&tt.stim.num==2),2),100))
    plot(MySmooth(mean(ts(:,tt.R&tt.hit&tt.stim.num==3),2),100))
    plot(MySmooth(mean(ts(:,tt.R&tt.hit&tt.stim.num==4),2),100))
end

figure(1);
for i=1:25
    unit = ix(i);
    ts = squeeze(psth(:,unit,:));
    subplot(5,5,i)
    hold on;
    plot(MySmooth(mean(ts(:,tt.L&tt.hit&tt.stim.num==0),2),100))
    plot(MySmooth(mean(ts(:,tt.L&tt.hit&tt.stim.num==1),2),100))
    plot(MySmooth(mean(ts(:,tt.L&tt.hit&tt.stim.num==2),2),100))
    plot(MySmooth(mean(ts(:,tt.L&tt.hit&tt.stim.num==3),2),100))
    plot(MySmooth(mean(ts(:,tt.L&tt.hit&tt.stim.num==4),2),100))
end




R = squeeze(mean(psth(:, :, tt.R&tt.hit&tt.stim.num==0), 3));
Rstim = squeeze(mean(psth(:, :, tt.R&tt.hit&tt.stim.num==1), 3));

Rdiff = Rstim-R;
figure; imagesc(MySmooth(Rdiff(:, ix(1:106)), 100)); caxis([-10 10]);
figure; plot(mean(MySmooth(Rdiff(:, ix(1:106)), 100), 2))




L = squeeze(mean(psth(:, :, tt.L&tt.hit&tt.stim.num==0), 3));
Lstim = squeeze(mean(psth(:, :, tt.L&tt.hit&tt.stim.num==1), 3));

Ldiff = Lstim-L;
figure; imagesc(MySmooth(Ldiff(:, ix(1:106)), 100)); caxis([-20 20]);
figure; plot(mean(MySmooth(Ldiff(:, ix(1:106)), 100), 2))





function cv = calcCodingVector(psth, meta, params)
% dt = 1./25000;
% time = (1:size(psth, 1))./dt;


Nclust = size(psth, 2);
mu = zeros(Nclust, 2);
sd = zeros(Nclust, 2);


for i = 1:Nclust
    ts1 = squeeze(psth(params.ix1, i, params.trial1));
    ts2 = squeeze(psth(params.ix2, i, params.trial2));
    
    mu(i, 1) = mean(mean(ts1));
    mu(i, 2) = mean(mean(ts2));
    
    sd(i, 1) = std(mean(ts1, 1), [], 2);
    sd(i, 2) = std(mean(ts2, 1), [], 2);
end
cv = (mu(:,2)-mu(:,1))./sqrt(sd(:,1).^2 + sd(:,2).^2);

probeix = ismember(meta.probe, params.probe);
cv(~probeix) = 0;

cv(isnan(cv)) = 0;
cv = cv./sum(abs(cv));




function plotCVProj(obj, psth, cv, trials, plotclr, trialsToPlot, sm)
% figure; hold on;
time = (1:size(psth, 1))./200;
proj = zeros(size(psth, 1), sum(trials(trialsToPlot)));
traj = zeros(10000, sum(trials(trialsToPlot)));
cnt = 0;


for i = trialsToPlot
    if ~trials(i)
        continue;
    end
    
    cnt = cnt+1;
    
    dat = MySmooth(psth(:, :, i), sm);
    proj(:, cnt) = dat*cv;
    traj(1:size(obj.traj{1}(i).ts, 1), cnt) = obj.traj{1}(i).ts(:, 2, 1);
end
% 
% plot(time, proj, 'Color', plotclr);
% plot(time, mean(proj, 2), 'Color', plotclr*0.5, 'LineWidth', 3);

figure; imagesc(proj)
figure; imagesc(traj(1:2000, :))










function plotPCAProj(psth, COEFF, trials, plotclr)

for i = 1:numel(trials)
    if trials(i)
        clr = plotclr;
    else
        continue;
    end
    
    dat = MySmooth(psth(:, :, i), 400);
    proj = dat*COEFF;
    
    ix = 10:70;
    plot3(proj(ix, 1), proj(ix, 2), proj(ix,3), 'Color', 'k', 'LineWidth', 2);
    
    ix = 70:520;
    plot3(proj(ix, 1), proj(ix, 2), proj(ix,3), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    ix = 520;
    plot3(proj(ix, 1), proj(ix, 2), proj(ix,3), '.', 'MarkerSize', 50, 'Color', [0 0 0], 'LineWidth', 2);
    
    plot3(proj(:, 1), proj(:, 2), proj(:,3), 'Color', clr);
    
end










function [tt, psth, meta, obj] = ObjToPSTH(obj, edges)

dt = mean(diff(edges));




tt = obj.bp;
tt.lickL = obj.bp.R&obj.bp.miss | obj.bp.L&obj.bp.hit;
tt.lickR = obj.bp.L&obj.bp.miss | obj.bp.R&obj.bp.hit;


tt.trialtypenames = unique(obj.bp.stim.num);
NtrialTypes = numel(tt.trialtypenames);
tt.trialtypeflags = cell(NtrialTypes, 1);

for i = 1:NtrialTypes
    
    tt.trialtypeflags{i} = obj.bp.stim.num==i;
    
end

Nunits = 0;
for i = 1:numel(obj.clu)
    for j = 1:numel(obj.clu{i})
        Nunits = Nunits+1;
    end
end



Ntrials = obj.bp.Ntrials;

psth = NaN.*zeros(numel(edges), Nunits, Ntrials);

quality = cell(Nunits, 1);
channel = zeros(Nunits, 1);
probe = zeros(Nunits, 1);

cnt = 0;
for i = 1:numel(obj.clu)
    NunitsOnProbe = numel(obj.clu{i});
    
    for j = 1:NunitsOnProbe
        cnt = cnt+1;
        
        unit = obj.clu{i}(j);
        
        spktrial = unit.trial;
        spktrialtm = unit.trialtm;

        trial1 = 1;
        trial2 = Ntrials;
        
        for k = trial1:trial2
            ix = spktrial==k;%obj.trialIDs(i);
            cuetm = obj.sglx.rewardFileOffset;
            t = spktrialtm(ix)-cuetm(k);
            
            psth(:, cnt, k) = histc(t, edges)./dt;
            
        end
        
        
        if isempty(obj.clu{i}(j).site)
            channel(cnt) = 1;
        else
            channel(cnt) = obj.clu{i}(j).site;
        end
        
        quality{cnt} = obj.clu{i}(j).quality;
        probe(cnt) = i;
    end
end
psth = psth(1:end-1, :, :);

meta.unitNumber = NunitsOnProbe;
meta.trialTypeStr = tt.trialtypenames;
% meta.trialOutcomeStr = obj.trialTypeStr;
% meta.trialOutcomeFlags = obj.trialTypeMat;
meta.unitQuality = quality;
meta.channel = channel;
% siteMap = getSiteLocs(obj.sessionMeta.probeType, obj.sessionMeta.depth);
% meta.depth = siteMap(meta.channel);
% meta.manipulatorDepth = obj.sessionMeta.depth;
meta.unitNum = Nunits;
meta.probe = probe;



function siteLocs = getSiteLocs(probeType, depth)


switch probeType
    case 'janelia2x32'
        siteLocs = [depth - 20 - (31:-1:0)*25 depth - 20 - (31:-1:0)*25]';        

end
















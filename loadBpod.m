function bp = loadBpod(pth, fn)

% pth = 'Y:\EEL\Experiments\EEL2\Bpod';
% fn = 'EEL2_soundTaskMultiLickBitcode_20190910_162350';

fn = fullfile(pth, fn);
dat = load(fn);

bp = unpackBpodData(dat);
% plotBpodData(bp);




function plotBpodData(bp)

figure; 
subplot(1,2,1); hold on;
plotByLick(bp, find(bp.L));

subplot(1,2,2); hold on;
plotByLick(bp, find(bp.R));





function plotByLick(bp, alltrials)
% figure; hold on;
stimtype = unique(bp.stim.num);
Nplotted = 0;
for i = 1:numel(stimtype)
    trials = alltrials(bp.stim.num(alltrials)==stimtype(i));
    
    y = [Nplotted Nplotted+numel(trials)];
    
    if stimtype(i)==1
        myfill([1.6 2.0], y, 0.2);
    elseif stimtype(i)==2
        myfill([1.9 2.3], y, 0.2);
    elseif stimtype(i)==3
        myfill([2.5 3.5], y, 0.2);
    elseif stimtype(i)==4
        myfill([3 4], y, 0.2);
    end
    
    for j = 1:numel(trials)
        
        lickL = bp.ev.lickL{trials(j)};
        lickR = bp.ev.lickR{trials(j)};
        
        plot(lickL, Nplotted+j*ones(size(lickL)), 'r.');
        plot(lickR, Nplotted+j*ones(size(lickR)), 'b.');
        
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
xlim([0 10]);
ylim([0 Nplotted]);


function myfill(xl, yl, rampdur)
x = [xl(1) xl(1) xl(2) xl(2)];
y = [yl(1) yl(2) yl(2) yl(1)];
f = fill(x, y, [0.7 0.7 1]); set(f, 'LineStyle', 'None');
f = fill([x(3:4) x(3:4)+rampdur], y, [0.85 0.85 1]); set(f, 'LineStyle', 'None');


function bp = unpackBpodData(dat)
dat = dat.SessionData;

bp.Ntrials = dat.nTrials;

bp.hit = NaN.*zeros(bp.Ntrials, 1);
bp.R = (dat.TrialTypes == 0)';
bp.L = (dat.TrialTypes == 1)';
bp.miss = NaN.*zeros(bp.Ntrials, 1);
bp.no = NaN.*zeros(bp.Ntrials, 1);
bp.early = NaN.*zeros(bp.Ntrials, 1);

bp.protocol.types = cell(bp.Ntrials, 1);
bp.protocol.nums = zeros(bp.Ntrials, 1);

bp.ev.bitStart = NaN.*zeros(bp.Ntrials, 1);
bp.ev.sample = NaN.*zeros(bp.Ntrials, 1);
bp.ev.delay = NaN.*zeros(bp.Ntrials, 1);
bp.ev.goCue = NaN.*zeros(bp.Ntrials, 1);
bp.ev.reward = NaN.*zeros(bp.Ntrials, 1);

bp.ev.lickL = cell(bp.Ntrials, 1);
bp.ev.lickR = cell(bp.Ntrials, 1);

bp.stim.enable = NaN.*zeros(bp.Ntrials, 1);
bp.stim.num = NaN.*zeros(bp.Ntrials, 1);
bp.stim.amp = cell(bp.Ntrials, 1);
bp.stim.del = cell(bp.Ntrials, 1);
bp.stim.dur = cell(bp.Ntrials, 1);
bp.stim.loc = cell(bp.Ntrials, 1);
bp.stim.state = cell(bp.Ntrials, 1);

% if iscell(dat.TrialParams)
%     if ~isempty(strfind(dat.TrialSettings(1).GUIMeta.Location.String{dat.TrialSettings(1).GUI.Location}, 'EphysRig1'))
%         bp.stim.amp{1} = dat.TrialParams{1}.wav.stim.amp;
%         bp.stim.state{1} = dat.TrialParams{1}.wav.stim.state;
%         bp.stim.del{1} = dat.TrialParams{1}.wav.stim.del;
%         bp.stim.dur{1} = dat.TrialParams{1}.wav.stim.dur;
%     end
% else
%     bp.stim.amp = NaN;
%     bp.stim.state = NaN;
%     bp.stim.del = NaN;
%     bp.stim.dur = NaN;
% end

bp.bitRand = NaN.*zeros(bp.Ntrials, 1);
bp.protocol.types = dat.TrialSettings(1).GUIMeta.ProtocolType.String;
bp.autowater.types = dat.TrialSettings(1).GUIMeta.Autowater.String;
bp.autolearn.types = dat.TrialSettings(1).GUIMeta.Autolearn.String;

for i = 1:bp.Ntrials
    

    
    if isfield(dat.RawEvents.Trial{i}.States, 'Reward')
        bp.hit(i) = ~isnan(dat.RawEvents.Trial{i}.States.Reward(1));
    else
        bp.hit(i) = NaN;
    end
    if isfield(dat.RawEvents.Trial{i}.States, 'TimeOut')
        bp.miss(i) = ~isnan(dat.RawEvents.Trial{i}.States.TimeOut(1));
    else
        bp.miss(i) = NaN;
    end
    if isfield(dat.RawEvents.Trial{i}.States, 'NoResponse')
        bp.no(i) = ~isnan(dat.RawEvents.Trial{i}.States.NoResponse(1));
    else
        bp.no(i) = NaN;
    end
    
    if isfield(dat.RawEvents.Trial{i}.States, 'EarlyLickSample')
        if isfield(dat.RawEvents.Trial{i}.States, 'EarlyLickDelay')
            bp.early(i) = ~isnan(dat.RawEvents.Trial{i}.States.EarlyLickSample(1)) | ~isnan(dat.RawEvents.Trial{i}.States.EarlyLickDelay(1));
        else
            bp.early(i) = ~isnan(dat.RawEvents.Trial{i}.States.EarlyLickSample(1));
        end
    end
    
    if isfield(dat.RawEvents.Trial{i}.States, 'bit1')
        bp.ev.bitStart(i) = dat.RawEvents.Trial{i}.States.bit1(1);
    end

    if isfield(dat.RawEvents.Trial{i}.States, 'SamplePeriod')
        bp.ev.sample(i) = dat.RawEvents.Trial{i}.States.SamplePeriod(1);
    end
    
    if isfield(dat.RawEvents.Trial{i}.States, 'DelayPeriod')
        bp.ev.delay(i) = dat.RawEvents.Trial{i}.States.DelayPeriod(1);
    end
    
    if isfield(dat.RawEvents.Trial{i}.States, 'GoCue')
        bp.ev.goCue(i) = dat.RawEvents.Trial{i}.States.GoCue(1);
    end
    
    if isfield(dat.RawEvents.Trial{i}.States, 'Reward')
        bp.ev.reward(i) = dat.RawEvents.Trial{i}.States.Reward(1);
    end
    
    if isfield(dat.RawEvents.Trial{i}.Events, 'Port1In')
        bp.ev.lickL{i} = dat.RawEvents.Trial{i}.Events.Port1In;
    end
    
    if isfield(dat.RawEvents.Trial{i}.Events, 'Port2In')
        bp.ev.lickR{i} = dat.RawEvents.Trial{i}.Events.Port2In;
    end
    
    if iscell(dat.TrialParams)
        bp.bitRand(i) = dat.TrialParams{i}.bit.randNum;%dat.TrialParams(i).bit.randNum;
    else
        bp.bitRand(i) = dat.TrialParams(i).bit.randNum;
    end
    
    bp.protocol.nums(i) = dat.TrialSettings(i).GUI.ProtocolType;
    bp.autolearn.nums(i) = dat.TrialSettings(i).GUI.Autolearn;
    bp.autowater.nums(i) = dat.TrialSettings(i).GUI.Autowater;
    
    if iscell(dat.TrialParams)
        
        bp.stim.enable(i) = dat.TrialParams{i}.giveStim;
        bp.stim.num(i) = dat.TrialParams{i}.stimNum;
        if bp.stim.enable(i)
            if isfield(dat.TrialParams{i}, 'wav')
                bp.stim.amp{i} = dat.TrialParams{i}.wav.stim.amp{bp.stim.num(i)};
            else
                bp.stim.amp{i} = NaN;
            end
            bp.stim.del{i} = dat.TrialParams{i}.stimDel;
%             bp.stim.dur{i} = dat.TrialParams{i}.stimDur;
%             bp.stim.loc{i} = dat.TrialParams{i}.stimLoc;
            bp.stim.state{i} = dat.TrialParams{i}.stimState;
        end
%         if strcmp(dat.Info.SessionDate, '27-Feb-2020')
%             if bp.stim.num(i) == 1
%                 bp.stim.num(i) = 3;
%             end
%         end
    else
        bp.stim.enable(i) = NaN;
        bp.stim.num(i) = NaN;
    end
 
end

% stimNums = unique(bp.stim.num); %stimNums is sorted
% for i=1:length(stimNums)
%     bp.stim.numParams.num = stimNums(i);
%     if stimNums(i)==0
%         bp.stim.numParams.amp = 0;
%         bp.stim.numParams.state = '';
%         bp.stim.numParams.del = 0;
%         bp.stim.numParams.dur = 0;
%     else
%         bp.stim.numParams.amp = bp.stim.amp(stimNums(i));
%         bp.stim.numParams.state = bp.stim.state{stimNums(i)};
%         bp.stim.numParams.del = bp.stim.del(stimNums(i));
%         bp.stim.numParams.dur = bp.stim.dur(stimNums(i));
%     end
% end




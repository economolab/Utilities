function createDataObject(pth)
obj.pth = pth;

obj.bp = getBpodData(pth);
obj.sglx = getSglxData(pth);
jrc = getJRClustData(pth);
obj.traj = getTrajectories(pth);

obj.trials = getTrials(pth, obj);

obj.clu = makeClu(obj, jrc);
obj = syncTrials(obj);    % TODO: better interface for manually correcting trials
obj = syncTimeBase(obj);  % TODO: which frames does the camera drop

save(fullfile(pth.sv, pth.fn.obj), 'obj');







function obj = syncTimeBase(obj)
%Notes obj.bp.ev:  t=0 start of 'TrialStart' state. Exampe timing  
%   bitcode start   :   t = 0.04
%   sample          :   t = 0.3
%   delay           :   t = 1.6
%   goCue           :   t = 2.5


%Align everything to start of bitcode
obj.bp.ev = subtractBitStart(obj.bp.ev, obj.bp.ev.bitStart);
obj.bp.ev.timeBase  = 'Seconds after BitCode start';

%Spike GLX timing
%   obj.sglx.bitFileOffset        :    ~0.53 (padded by 0.5 sec. seems like should be 0.51) 
%   obj.sglx.cueFileOffset        :     2.99 (formerly rewardFileOffset)

obj.sglx.bitFileOffset_uncorrected = obj.sglx.bitFileOffset;
for i = 1:numel(obj.sglx.meta.camTrigIX)
    obj.sglx.meta.camTrigIX{i} = obj.sglx.meta.camTrigIX{i} - (obj.sglx.bitFileOffset(i)*obj.sglx.meta.fs);
end
obj.sglx.cueFileOffset = obj.sglx.cueFileOffset - obj.sglx.bitFileOffset;
obj.sglx.bitFileOffset = obj.sglx.bitFileOffset - obj.sglx.bitFileOffset;
obj.sglx.timeBase = 'Seconds after BitCode start';


%Align Neural data to start of bitcode.
for j = 1:numel(obj.clu)
    for i = 1:numel(obj.clu{j})
        obj.clu{j}(i).trialtm = obj.clu{j}(i).tm - obj.sglx.fileStart(obj.clu{j}(i).trial);
        obj.clu{j}(i).trialtm = obj.clu{j}(i).trialtm - obj.sglx.bitFileOffset(i);
        obj.clu{j}(i).timeBase = 'Seconds after BitCode start';
    end
end


%Camera timing
%   Frame number 1 should be triggered 20 ms before start of bitcode.  An
%   indeterminate number of frames are dropped at the (presumably at the 
%   beginning; usually 2).  Number of dropped frames can be calculated.


%Assuming dropped frames all at the beginning
FrameRate = 400;
Ntrigs = cellfun(@numel, obj.sglx.meta.camTrigIX);

for i = 1:numel(obj.traj)
    for j = 1:numel(obj.traj{i})
        df = Ntrigs(j) - size(obj.traj{i}(j).ts, 1);
        obj.traj{i}(j).time = (df+((0:size(obj.traj{i}(j).ts, 1)-1)))./FrameRate - 0.02; %20 ms before bitcode start
        obj.traj{i}(j).timeBase = 'Seconds after BitCode start';
        obj.traj{i}(j).droppedFrames = df;
    end
end











function obj = syncTrials(obj)

% modify bp stuff
usable = getUsableTrials(obj);

%Bpod
bp = obj.bp;
bp = subsetStruct(bp, usable.bpod, bp.Ntrials);
bp.Ntrials = sum(usable.bpod);
obj.bp = bp;

%SGLX
sglx = obj.sglx;
sglx = subsetStruct(sglx, usable.sglx(usable.sglx>0), numel(sglx.fileStart));
obj.sglx = sglx;

%Video
traj = obj.traj;
for i = 1:numel(traj)
    a.traj = traj{i};
    a = subsetStruct(a, usable.vid(usable.vid>0), numel(a.traj));
    traj{i} = a.traj;
end
obj.traj = traj;

%Trials
trials = obj.trials;
trials = subsetStruct(trials, usable.trials, obj.trials.N);
obj.trials = trials;


%Clusters
for i = 1:numel(obj.clu) % probes
    for j = 1:numel(obj.clu{i}) % units
        unit = obj.clu{i}(j);
        
        [foundspks, newTrials] = ismember(unit.trial, find(usable.sglx));
        unit.trial = newTrials;
        unit.trial(~foundspks) = NaN;
        
%         unit.trial = unit.trial(foundspks);
        
        obj.clu{i}(j) = unit;
    end
end




function s = subtractBitStart(s, val)

fnames = fieldnames(s);

N = numel(val);

for i= 1:numel(fnames)
    if numel(s.(fnames{i})) == N && isnumeric(s.(fnames{i}))
        s.(fnames{i}) = s.(fnames{i}) - val;
    elseif numel(s.(fnames{i})) == N && iscell(s.(fnames{i}))
        for j = 1:N
            if ~isempty(s.(fnames{i})(j))
                s.(fnames{i}){j} = s.(fnames{i}){j} - val(j);
            end
        end
    end
end



function s = subsetStruct(s, usable, N)

fnames = fieldnames(s);

for i= 1:numel(fnames)
    if numel(s.(fnames{i})) == N
        dat = s.(fnames{i});
        s.(fnames{i}) = dat(usable);
    elseif numel(s.(fnames{i}))==1 && isstruct(s.(fnames{i}))
        s.(fnames{i}) = subsetStruct(s.(fnames{i}), usable, N);
    end
end









function usable = getUsableTrials(obj)

usable.bpod = false(obj.trials.N, 1);
usable.trials = false(obj.trials.N, 1);
usable.sglx = zeros(obj.trials.N, 1);
usable.vid = zeros(obj.trials.N, 1);

cnt = 0;
for i = 1:obj.trials.N
    btrial = i;
    
    %do we have a SGLX file?
    if obj.trials.sglxNum(btrial)>0
        strial = obj.trials.sglxNum(btrial);
        if obj.trials.SGLXTrialsWithVid(strial)
            cnt = cnt+1;
            vtrial = obj.trials.VidFilesToUse(cnt);
            usable.bpod(btrial) = true;
            usable.sglx(btrial) = strial;
            usable.vid(btrial) = vtrial;
            usable.trials(btrial) = true;

        end
    end
    
end



function trials = getTrials(pth, obj)
trials.N = obj.bp.Ntrials;

if ~exist(fullfile(pth.sv, pth.fn.sglx), 'file')
    disp(['Could not find SpikeGLX file!  Expected to find: ' fullfile(pth.sv, pth.fn.sglx)]);
    return;
end
if ~exist(fullfile(pth.sv, pth.fn.bpod), 'file')
    disp(['Could not find Bpod file!  Expected to find: ' fullfile(pth.sv, pth.fn.bpod)]);
    return;
end

[trials.haveEphys, trials.sglxNum] = ismember((1:trials.N)', obj.sglx.meta.bitcode.trial);
[trials.haveBpod, trials.bpodNum] = ismember(obj.sglx.meta.bitcode.trial, (1:trials.N)');

Nviews = numel(obj.traj);
for i = 1:Nviews
    Nmovies = numel(obj.traj{i});
    for j = 1:Nmovies
        Nframes(j, i) = size(obj.traj{i}(j).ts, 1);
    end
end

if ~sum(diff(Nframes, 1, 2))
    disp('Frame check:  Cameras reporting same number of frames!');
else
    figure; plot(diff(Nframes, 1, 2)); title('Difference in frames between cameras');
    disp('!!!!Variable frame numbers across cameras');
    
end

vidFrames = mean(Nframes, 2);

Ntrigs = cellfun(@numel, obj.sglx.meta.camTrigIX);

C = xcorr(vidFrames-mean(vidFrames),Ntrigs-mean(Ntrigs));
[~, shift] = max(C);
shift = numel(Ntrigs)-shift;

SGLXTrialsToUse = 1:numel(Ntrigs);
VidTrialsToUse = 1:numel(vidFrames);
if ~shift
    disp('Have video for every spikeGLX file!')
elseif numel(vidFrames)~=Ntrigs
    DONE = 0;
    
    while ~DONE
       
        newTrials = inputdlg({['SGLXTrialsToUse (' num2str(numel(Ntrigs)) ')'], ['VidTrialsToUse (' num2str(numel(vidFrames)) ')']}...
            , 'Refine', 1, {['1:' num2str(numel(Ntrigs))], ['1:' num2str(numel(vidFrames))]});
        
        SGLXTrialsToUse = str2num(newTrials{1});
        VidTrialsToUse = str2num(newTrials{2});
        
        if numel(SGLXTrialsToUse)==numel(VidTrialsToUse)
            
            fig = figure;  hold on;
            plot(Ntrigs(SGLXTrialsToUse), 'r.-');
            plot(vidFrames(VidTrialsToUse), 'b');
            
            legend('SpikeGLX triggers', 'Video file frames');
            title(['Possible shift of ' num2str(shift) 'trials']);
            
            answer = questdlg('Accept/Refine results', 'Check', 'Accept', 'Refine', 'Accept');
            close(fig);
            
            if strcmp(answer, 'Accept')
                DONE=1;
            end
        end
        
    end
end

trials.SGLXTrialsWithVid = false(numel(Ntrigs), 1);
trials.SGLXTrialsWithVid(SGLXTrialsToUse) = true;

trials.VidFilesToUse = VidTrialsToUse';






function clu = makeClu(obj, jrc)

if isempty(jrc)
    clu = [];
    return
end

for j = 1:numel(jrc)
    Nclust = numel(jrc{j}.unitCount);
    % clu = zeros(Nclust, 1); %this might not work
    disp(['Adding ' num2str(Nclust) ' units from probe #' num2str(j)]);
    excellent = 0;
    good = 0;
    for i = 1:Nclust
        clu{j}(i).ix = jrc{j}.spikeTimes(jrc{j}.spikesByCluster{i});
        clu{j}(i).tm = double(clu{j}(i).ix)./obj.sglx.meta.fs;
        [~,~,sglxtrial] = histcounts(clu{j}(i).tm,[obj.sglx.fileStart' inf]);
        
        badspikes = obj.trials.bpodNum(sglxtrial)==0;
        clu{j}(i).ix = clu{j}(i).ix(~badspikes);
        clu{j}(i).tm = clu{j}(i).tm(~badspikes);
        sglxtrial = sglxtrial(~badspikes);
        
        clu{j}(i).trial = sglxtrial;
%         clu{j}(i).trialtm = clu{j}(i).tm - obj.sglx.fileStart(clu{j}(i).trial);
%         clu{j}(i).trialtm = clu{j}(i).trialtm - (obj.sglx.bitFileOffset(i) - obj.bp.ev.bitStart(obj.trials.bpodNum(i)));
        
        clu{j}(i).trial = obj.trials.bpodNum(clu{j}(i).trial);
        
        clu{j}(i).site = mode(jrc{j}.spikeSites(jrc{j}.spikesByCluster{i}));
        clu{j}(i).quality = jrc{j}.clusterNotes{i};
        if strcmp(clu{j}(i).quality, 'Excellent')
            excellent = excellent+1;
        elseif strcmp(clu{j}(i).quality, 'Good')
            good = good+1;
        end
    end
    disp(['      ' num2str(excellent) ' excellent units']);
    disp(['      ' num2str(good) ' good units']);
end


function  jrc = getJRClustData(pth)
contents = dir(pth.sglx);
good = false(numel(contents), 1);
for i = 1:numel(contents)
    if ~isempty(strfind(contents(i).name, '_res.mat'))
        good(i) = true;
        
    end
end
contents = contents(good);

for i = 1:numel(contents)
    jrc{i} = load(fullfile(pth.sglx, contents(i).name));
end

if ~exist('jrc', 'var')
    jrc = [];
end
   




function bp = getBpodData(pth)
fn.bp = fullfile(pth.sv, pth.fn.bpod);
if exist(fn.bp, 'file')
    a = load(fn.bp);
else
    disp(['WARNING: Could not find: ' fn.bp]);
    a.bp = [];
end
bp = a.bp;





function traj = getTrajectories(pth)
fn.vid = fullfile(pth.sv, pth.fn.vid);
if exist(fn.vid, 'file')
    a = load(fn.vid);
else
    disp(['WARNING: Could not find: ' fn.vid]);
    a.traj = [];
end
traj = a.traj;



function sglx = getSglxData(pth)
fn.sglx = fullfile(pth.sv, pth.fn.sglx);
if exist(fn.sglx, 'file')
    load(fn.sglx);
    
    sglx.fileStart = [0; cumsum(meta.Nsamp(1:end-1))./meta.fs];
    sglx.bitFileOffset = meta.bitcode.bitstart/25000;
    sglx.cueFileOffset = meta.bitcode.bitend/25000;
    sglx.meta = meta;
else
    disp(['WARNING: Could not find: ' fn.sglx]);
    sglx = [];
end



function plotTrialAverage(obj, clu, trialix, clr)
trialnum = find(trialix);
spks = ismember(obj.clu(clu).trial, trialnum);
spktms = obj.clu(clu).trialtm(spks);

edges = 0:0.05:5;
N = histc(spktms, edges);

binwid = mean(diff(edges));
plot(edges+binwid/2, N./numel(trialnum)./binwid, clr);








































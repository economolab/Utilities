function parseBpodData(pth, dt)

fns = findBpodFiles(pth.bpod, dt);

for i = 1:numel(fns)
    bp(i) = loadBpod(pth.bpod, fns{i});
end

if numel(bp)>1
    disp(['!!!Found ' num2str(numel(bp)) ' matching bpod files.  Concatenating them']);
    bp = concatenateBPfiles(bp);

end
save(fullfile(pth.sv, 'bpod.mat'), 'bp');




function bpall = concatenateBPfiles(bp)
bpall = bp(1);
for i = 2:numel(bp)
    bpall.Ntrials = bpall.Ntrials+bp(i).Ntrials;
    
    fnames = fieldnames(bp(i));
    for j = 1:numel(fnames)
        if ~strcmp(fnames{j}, 'ev') && ~strcmp(fnames{j}, 'Ntrials') && ~strcmp(fnames{j}, 'protocol') && ~strcmp(fnames{j}, 'stim')
            bpall.(fnames{j})(end+1:end+bp(i).Ntrials) = bp(i).(fnames{j});
        end
    end
    
end

for i = 2:numel(bp)
    
    fnames = fieldnames(bp(i).ev);
    for j = 1:numel(fnames)
        bpall.ev.(fnames{j})(end+1:end+bp(i).Ntrials) = bp(i).ev.(fnames{j});
    end
    fnames = fieldnames(bp(i).protocol);
    for j = 1:numel(fnames)
        if ~iscell(eval(['bp(i).protocol.' fnames{j}]))
            bpall.protocol.(fnames{j})(end+1:end+bp(i).Ntrials) = bp(i).protocol.(fnames{j});
        end
    end
    fnames = fieldnames(bp(i).stim);
    for j = 1:numel(fnames)
        if ~iscell(eval(['bp(i).stim.' fnames{j}]))
            bpall.stim.(fnames{j})(end+1:end+bp(i).Ntrials) = bp(i).stim.(fnames{j});
        end
    end
end
    




function fns = findBpodFiles(pth, dt)
contents = dir(pth);
good = false(size(contents));

for i = 1:numel(contents)
    day = dt(dt~='-');
    if contains(contents(i).name, day)
        good(i) = true;
    end
end
contents = contents(good);
fns = cell(numel(contents), 1);
for i = 1:numel(fns)
    fns{i} = contents(i).name;
end

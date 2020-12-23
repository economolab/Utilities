function skimSpikeGLXBinFiles(pth, svpth)

meta.pth = pth;
meta.fns = getMetaFiles(pth);

meta.fs = getMetaData(fullfile(meta.pth, meta.fns{1}), 'niSampRate');
meta.padSec = getMetaData(fullfile(meta.pth, meta.fns{1}), 'trgTTLMarginS');
meta.Nchan = getMetaData(fullfile(meta.pth, meta.fns{1}), 'nSavedChans');
meta.bytesPerSample = 2;

meta.bitcodeChan = 71;
meta.camTrigChan = 65;

meta = getBinMetaData(meta);

events = getGoCue(meta);
save(fullfile(svpth, 'events.mat'), 'events');
save(fullfile(svpth, 'MetaData.mat'), 'meta');




function meta = getBinMetaData(meta)
Nfiles = numel(meta.fns);
chans = [meta.bitcodeChan meta.camTrigChan];

meta.bitcode = cell(Nfiles, 1);
meta.camTrigIX = cell(Nfiles, 1);
meta.Nsamp = zeros(Nfiles, 1);

for i = 1:Nfiles
    disp(['(' num2str(i) '/' num2str(Nfiles) ')  Parsing file ' meta.fns{i}]);
    data = getBinData(meta, chans, i);
    meta.Nsamp(i) = size(data, 1);
    meta.bitcode{i} = parseBitCode(meta, data(:,1));
    meta.camTrigIX{i} = parseCamTrig(data(:,2));
end





function events = getGoCue(meta)
Ntrials = numel(meta.fns);
events = NaN.*zeros(Ntrials, 1);

events(1) = (meta.bitcode{1}.highIX(end))./meta.fs;
for i = 1:Ntrials-1
    if ~isempty(meta.bitcode{i}.highIX)
        events(i+1) = (sum(meta.Nsamp(1:i)) + meta.bitcode{i}.highIX(end))./meta.fs;
    end
end
events = events(~isnan(events));






function trig = parseCamTrig(data)
data = abs(data);
thresh = 2000;
trig = find(data(1:end-1)<thresh & data(2:end)>thresh)+1;




function code = parseBitCode(meta, data)
data = abs(data);

thresh = 2000;
bitRand = 2:13;
bitTrial = 14:25;
Nbits = numel(bitRand)+numel(bitTrial)+1;

code.highIX = [];
code.bits= zeros(1, Nbits);
code.rand = -1;
code.trial = -1;

try
    code.highIX = find(data(1:end-1)<thresh & data(2:end)>thresh)+1;
    bitlen = meta.fs*0.01;
    
    bitdat = reshape(data(code.highIX(1):code.highIX(1)+bitlen*Nbits-1), bitlen, Nbits);
    
    code.bits = mean(bitdat, 1)>thresh/2;
    code.rand = bin2dec(num2str(code.bits(bitRand)));
    code.trial = bin2dec(num2str(code.bits(bitTrial)));
catch
    disp('WARNING: Could not parse Bit Code!')
end





function data = getBinData(meta, chans, filenum)

Nsamp = getMetaData(fullfile(meta.pth, meta.fns{filenum}), 'fileSizeBytes')/meta.bytesPerSample/meta.Nchan;

data = zeros(Nsamp, numel(chans));

for i = 1:numel(chans)
    fn = fullfile(meta.pth, [meta.fns{filenum}(1:end-5) '.bin']);
    fid = fopen(fn, 'r');
    fseek(fid, meta.bytesPerSample*(chans(i)-1), 0);
    data(:,i) = fread(fid, Nsamp, 'int16', meta.bytesPerSample*(meta.Nchan-1));
    fclose(fid);
end






function fn = getMetaFiles(pth)
contents = dir(pth);

good = false(size(contents));

for i = 1:numel(contents)
    if ~isempty(strfind(contents(i).name, '.meta'))
        good(i) = true;
    end
end

contents = contents(good);
datenum =  [contents(:).datenum];
[~, ix] = sort(datenum, 'ascend');
contents = contents(ix);

fn = cell(numel(contents), 1);
for i = 1:numel(contents)
    fn{i} = contents(i).name;
end

    


function val = getMetaData(fn, token)

fid = fopen(fn);

while ~feof(fid)
    
    str = fgetl(fid);
    if contains(str, token)
        pos = strfind(str, '=');
        val = str2double(str(pos+1:end));
    end
    
end
fclose(fid);




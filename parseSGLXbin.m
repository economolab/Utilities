function parseSGLXbin(pth)

meta.pth = pth.sglx;
meta.fns = getMetaFiles(pth.sglx);



meta.fs = getMetaData(fullfile(meta.pth, meta.fns{1}), 'niSampRate');
meta.padSec = getMetaData(fullfile(meta.pth, meta.fns{1}), 'trgTTLMarginS');
meta.Nchan = getMetaData(fullfile(meta.pth, meta.fns{1}), 'nSavedChans');

meta.bytesPerSample = 2;

meta.bitcodeChan = meta.Nchan-31;  %161 for 10/22, 97 for 10/23
meta.camTrigChan = meta.Nchan-30;  %162 for 10/22, 98 for 10/23


concatFn = fullfile(meta.pth, 'ConcatenatedData_Probe1.bin');
if exist(concatFn, 'file')
    answer = questdlg(['Delete ' concatFn]);
    if strcmp(answer, 'Yes')
        delete(concatFn);
    end
end
concatFn = fullfile(meta.pth, 'ConcatenatedData_Probe2.bin');
if exist(concatFn, 'file')
    answer = questdlg(['Delete ' concatFn]);
    if strcmp(answer, 'Yes')
        delete(concatFn);
    end
end


meta = getBinMetaData(meta);

events = getGoCue(meta);
save(fullfile(pth.sv, 'events.mat'), 'events');
save(fullfile(pth.sv, pth.fn.sglx), 'meta');




function meta = getBinMetaData(meta)
Nfiles = numel(meta.fns);
chans = [meta.bitcodeChan meta.camTrigChan];

val.bitcode = cell(Nfiles, 1);
val.camTrigIX = cell(Nfiles, 1);
val.Nsamp = zeros(Nfiles, 1);
val.noise = zeros(Nfiles, meta.Nchan);

Nprobes = (meta.Nchan==128) + 2*(meta.Nchan==192);

concatenate = 1;
calcNoise = 1;
for i = 1:Nfiles
    disp(['(' num2str(i) '/' num2str(Nfiles) ')  Parsing file ' meta.fns{i}]);
    data = getBinData(meta, chans, i); %0.38uV/bin
    val.Nsamp(i) = size(data, 1);
    val.bitcode{i} = parseBitCode(meta, data(:,1));
    val.camTrigIX{i} = parseCamTrig(data(:,2));
    
    disp(['Read bitcode ' num2str(val.bitcode{i}.trial) ]);

    if calcNoise
        siteMap = [47, 45, 44, 42, 40, 38, 36, 31, 33, 35, 37, 39, 41, 43, 46, 48, 21, 23, 25, 30, 28, 17, 19, 32, 20, 18, 34, 27, 29, 26, 24, 22, 12, 10, 8, 3, 5, 16, 14, 1, 13, 15, 63, 6, 4, 7, 9, 11, 50, 52, 53, 55, 57, 59, 61, 2, 64, 62, 60, 58, 56, 54, 51, 49];
        data = getAllBinData(meta, siteMap, i, Inf);
        [B,A] = butter(2,0.0240, 'high');
        data = filtfilt(B, A, data);
        val.noise(i, 1:size(data, 2)) = std(data, [], 1);
    end
    
    if concatenate
        siteMap = [47, 45, 44, 42, 40, 38, 36, 31, 33, 35, 37, 39, 41, 43, 46, 48, 21, 23, 25, 30, 28, 17, 19, 32, 20, 18, 34, 27, 29, 26, 24, 22, 12, 10, 8, 3, 5, 16, 14, 1, 13, 15, 63, 6, 4, 7, 9, 11, 50, 52, 53, 55, 57, 59, 61, 2, 64, 62, 60, 58, 56, 54, 51, 49];
        if Nprobes==2
            siteMap = [siteMap siteMap+64];
        end
        dataAll = getAllBinData(meta, siteMap, i, Inf);
        for k = 2%1:Nprobes %probes
            
            [B,A] = butter(3,300./25000, 'high');
            data = filtfilt(B, A, dataAll(:, 64*(k-1)+1:64*k));
            
%             for j = 1:size(data, 2)
%                 figure(floor((j-1)/32)+1); subplot(4, 8, j-32*floor((j-1)/32));
%                 plot(data(:, j)); ylim([-1500 500]); hold on; plot(mean(data, 2));
%             end


            dataNew = data;
%             inds = {1:32};
%             for jj = 1
%                 for j = 1:numel(inds{jj})
%                     
%                     use = (j-4:j+4);
%                     use = use(use>0&use<=32&use~=j);
%                     use = inds{jj}(use);
%                     dataNew(:, inds{jj}(j)) = data(:, inds{jj}(j)) - mean(data(:, use), 2);
%                 end
%             end
            dataNew = dataNew - median(dataNew, 2);
            
            figure; hold on;
            time = (1:size(dataNew, 1))/25000;
            
            for j = 1:5%size(dataNew, 2)
                plot(time, data(:, j));%+j*400);
            end

            for j = 1:5%size(dataNew, 2)
                plot(time, dataNew(:, j)+j*400);
            end
            
            data = dataNew; %0.38uV/bin
            
%             figure; plot(data(:, 1:32), 'k');
%             hold on; plot(data(:, 2:2:32), 'b');
%              hold on; plot(mean(data(:, 1:2:32), 2), 'r');
%             data(:, 1:32) = data(:, 1:32) - repmat(mean(data(:, 1:32), 2), 1, size(data(:, 1:32), 2));
%             data(:, 2:2:32) = data(:, 2:2:32) - repmat(mean(data(:, 2:2:32), 2), 1, size(data(:, 2:2:32), 2));
%             hold on; plot(data(:, 1:32), 'g');


            


% hold on; plot(dataNew(:, 1:32), 'r');



            thresh = 1e5;
            pow = var(data, [], 2);
%             pow = mean(dataNew.^2, 2); %Maybe try this instead
            ix = imdilate(pow>thresh, ones(13, 1));
            data(ix, :) = 0;
            

            
%             appendToBinaryFile(meta, data, k);
        end
        
    end
end

meta.bitcode.trial = zeros(Nfiles, 1);
meta.bitcode.rand = zeros(Nfiles, 1);
meta.bitcode.bitstart = zeros(Nfiles, 1);
meta.bitcode.bitend = zeros(Nfiles, 1);
meta.bitcode.bits = cell(Nfiles, 1);
meta.bitcode.highIX = cell(Nfiles, 1);


for i = 1:Nfiles
    meta.bitcode.trial(i) = val.bitcode{i}.trial;
    meta.bitcode.rand(i) = val.bitcode{i}.rand;
    
    meta.bitcode.bits{i} = val.bitcode{i}.bits;
    meta.bitcode.highIX{i} = val.bitcode{i}.highIX;
    
    if ~isempty(meta.bitcode.highIX{i})
        meta.bitcode.bitstart(i) = meta.bitcode.highIX{i}(1);
        meta.bitcode.bitend(i) = meta.bitcode.highIX{i}(end);
    end
end

meta.noise = val.noise;
meta.Nsamp = val.Nsamp;
meta.camTrigIX = val.camTrigIX;




function appendToBinaryFile(meta, data, num)

fid = fopen(fullfile(meta.pth, ['ConcatenatedData_Probe' num2str(num) '.bin']), 'a');
data = data';
fwrite(fid, data(:), 'int16');
fclose(fid);








function events = getGoCue(meta)
Ntrials = numel(meta.fns);
events = NaN.*zeros(Ntrials, 1);

events(1) = (meta.bitcode.highIX{1}(end))./meta.fs;
for i = 1:Ntrials-1
    if ~isempty(meta.bitcode.highIX{i})
        events(i+1) = (sum(meta.Nsamp(1:i)) + meta.bitcode.highIX{i}(end))./meta.fs;
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



function data = getAllBinData(meta, chans, filenum, npts)

fn = fullfile(meta.pth, [meta.fns{filenum}(1:end-5) '.bin']);
fid = fopen(fn, 'r');
data = fread(fid, meta.Nchan*npts, 'int16');
fclose(fid);
data = reshape(data, meta.Nchan, numel(data)./meta.Nchan)';
data = data(:, chans);





function data = getBinData(meta, chans, filenum)

Nsamp = getMetaData(fullfile(meta.pth, meta.fns{filenum}), 'fileSizeBytes')/meta.bytesPerSample/meta.Nchan;

data = zeros(Nsamp, numel(chans));

for i = 1:numel(chans)
    fn = fullfile(meta.pth, [meta.fns{filenum}(1:end-5) '.bin']);
    fid = fopen(fn, 'r');
    fseek(fid, meta.bytesPerSample*(chans(i)-1), 0);
    a = fread(fid, Nsamp, 'int16', meta.bytesPerSample*(meta.Nchan-1));
    data(1:numel(a),i) = a;
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




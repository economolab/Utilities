function parseSGLXbin_v2(pth, params)

meta.pth = pth.sglx;
meta.fns = getMetaFiles(pth.sglx);



meta.fs = getMetaData(fullfile(meta.pth, meta.fns{1}), 'niSampRate');
meta.padSec = getMetaData(fullfile(meta.pth, meta.fns{1}), 'trgTTLMarginS');
meta.Nchan = getMetaData(fullfile(meta.pth, meta.fns{1}), 'nSavedChans');

meta.bytesPerSample = 2;

meta.bitcodeChan = meta.Nchan-24;  %168 for 08/01, 161 for 10/22, 97 for 10/23
meta.camTrigChan = meta.Nchan-30;  %162 for 08/01, 162 for 10/22, 98 for 10/23


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


meta = parseBinData(meta, params);

events = getGoCue(meta);
save(fullfile(pth.sv, 'events.mat'), 'events');
save(fullfile(pth.sv, pth.fn.sglx), 'meta');




function meta = parseBinData(meta, params)
Nfiles = numel(meta.fns);
chans = [meta.bitcodeChan meta.camTrigChan];

val.bitcode = cell(Nfiles, 1);
val.camTrigIX = cell(Nfiles, 1);
val.Nsamp = zeros(Nfiles, 1);
val.noise = zeros(Nfiles, meta.Nchan);

val.Nprobes = (meta.Nchan==128) + 2*(meta.Nchan==192);

siteMap = params.siteMap;
for i = 2:val.Nprobes
    siteMap(end+1:end+64) = siteMap+64*(i-1);
end

for i = 1:Nfiles
    disp(['(' num2str(i) '/' num2str(Nfiles) ')  Parsing file ' meta.fns{i}]);
    data = getBinData(meta, chans, i); %0.38uV/bin
    val.Nsamp(i) = size(data, 1);
    val.bitcode{i} = parseBitCode(meta, data(:,1));
    val.camTrigIX{i} = parseCamTrig(data(:,2));
    
    disp(['Read bitcode ' num2str(val.bitcode{i}.trial) ]);
    
    if params.concatenate

        dataAll = getAllBinData(meta, siteMap, i, Inf);
        for k = 1:val.Nprobes
            
            [B,A] = butter(3,300./25000, 'high');
            data = filtfilt(B, A, dataAll(:, 64*(k-1)+1:64*k));

            dataMed = data;
            
            inds = {1:32; 33:64};
            for jj = 1:numel(inds)
                if params.medfilt>0
                    for j = 1:numel(inds{jj})
                        
                        use = (j-floor(params.medfilt/2):j+floor(params.medfilt/2));
                        use = use(use>0&use<=32&use~=j);
                        use = inds{jj}(use);
                        
                        dataMed(:, inds{jj}(j)) = data(:, inds{jj}(j)) - median(data(:, use), 2);
                    end
                else
                    dataMed(:, inds{jj}) = data(:, inds{jj}) - median(data(:, inds{jj}), 2);
                end
            end

            dataMedBlank = dataMed;
            
            pow = var(dataMedBlank, [], 2);
%             pow = mean(dataNew.^2, 2); %Maybe try this instead
            
            ix = imdilate(pow>params.powThresh, ones(13, 1));
            dataMedBlank(ix, :) = 0;
            
            if params.showPlots
                
                figure(100);
                hold off;
                plot(pow);
                hold on;
                plot([1 numel(pow)], params.powThresh+[0 0], 'k:');
                
                for j = 1:size(data, 2)
                    figure(floor((j-1)/32)+1); subplot(4, 8, j-32*floor((j-1)/32));  hold off; 
                    plot(data(:, j).*params.uVperBit); ylim([-500 250]); 
                    hold on; 
                    plot(dataMed(:,j).*params.uVperBit);
                end
                
                figure(200); 
                time = (1:size(dataMedBlank, 1))/25000;
                hold off;
                for j = 1:size(dataMedBlank, 2)
                    
                    plot(time, dataMedBlank(:, j)+j*400);
                    hold on;
                end
                axis tight;
                
                pause;
            end
            
            appendToBinaryFile(meta, dataMedBlank, k);
            
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




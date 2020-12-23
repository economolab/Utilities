function expPipeline()
% dataPath = 'Y:\EEL\Experiments\EEL2\';
% dates = {'2019-09-10'; '2019-09-11'; '2019-09-12'; '2019-09-13'; '2019-09-14'};

% dataPath = {'V:\EEL\Experiments\EEL1';}
% dates = {'2019-10-22', '2019-10-23', '2019-10-24'};
% dates = {'2019-10-25'};% };

% dataPath = 'V:\EEL\Experiments\EEL8';
% dates = {'2019-11-22'};
% % 
% dataPath = {'V:\EEL\Experiments\EEL11'};
% dates = {'2020-09-04'};
% dates = {'2020-03-01'};

% dataPath = {'E:\EEL7'};
% dates = {'2020-02-27','2020-02-28','2020-02-29','2020-03-01','2020-03-02'};
% dates = {'2020-03-01'};

% dataPath = {'Y:\EEL\Experiments\EEL14'};
%{'2020-08-01','2020-08-02','2020-08-03',}
% dates = {'2020-08-05','2020-08-07','2020-08-08',};
% dates = {'2020-08-09','2020-08-10','2020-08-11'};

% error in file t360 '2020-08-04'
% Index in position 2 exceeds array bounds.
% 
% Error in parseSGLXbin_v2>parseBinData (line 87)
%                         dataMed(:, inds{jj}(j)) = data(:, inds{jj}(j)) - median(data(:, use), 2);
% 
% Error in parseSGLXbin_v2 (line 34)
% meta = parseBinData(meta, params);
% 
% Error in expPipeline (line 37)
%         parseSGLXbin_v2(pth, SGLXparams);



% dataPath = {'Y:\EEL\Experiments\EEL15'};
% dates = {'2020-08-12', '2020-08-13', '2020-08-14', '2020-08-15', '2020-08-16', '2020-08-19', '2020-08-20', '2020-08-21'}; 

dataPath = {'Y:\EEL\Experiments\EEL16'};
dates = {'2020-08-06'};

% dataPath = {'Y:\JEB\Experiments\JEB1'};
% dates={'2020-09-09', '2020-09-10', '2020-09-11', '2020-09-12', '2020-09-14','2020-09-15'};


SGLXparams = getSGLXparams();

for j = 1:numel(dataPath)
    for i = 1:numel(dates)
        pth = getPths(dataPath{j}, dates{i});
        
%             loadVideoTrajectories(pth)
%             parseBpodData(pth, dates{i})
        
        parseSGLXbin_v2(pth, SGLXparams);
%         createDataObject(pth);
    end
end



function params = getSGLXparams()
params.siteMap = [47, 45, 44, 42, 40, 38, 36, 31, 33, 35, 37, 39, 41, 43, 46, 48, 21, 23, 25, 30, 28, 17, 19, 32, 20, 18, 34, 27, 29, 26, 24, 22, 12, 10, 8, 3, 5, 16, 14, 1, 13, 15, 63, 6, 4, 7, 9, 11, 50, 52, 53, 55, 57, 59, 61, 2, 64, 62, 60, 58, 56, 54, 51, 49];
params.concatenate = 1;
params.medfilt = 7; %channels to calculte median over (if 0, then global per shank)
params.powThresh = 1e5; %For artifact blanking
params.showPlots = 0;
params.uVperBit =  0.38;




function pth = getPths(parent, dt)
anm = strsplit(parent, filesep);
anm = anm(~cellfun(@isempty, anm));
anm = anm{end};

if ~exist(fullfile(parent, 'Analysis'), 'dir')
    mkdir(fullfile(parent, 'Analysis'));
end

if ~exist(fullfile(parent, 'Analysis', dt), 'dir')
    mkdir(fullfile(parent, 'Analysis', dt));
end

pth.sv = fullfile(parent, 'Analysis', dt);
pth.sglx = fullfile(parent, 'SpikeGLX', dt);
pth.bpod = fullfile(parent, 'Bpod');
pth.vid = {fullfile(parent, 'Video', dt, 'Cam0'); fullfile(parent, 'Video', dt, 'Cam1');};

pth.fn.sglx = 'SGLXmeta.mat';
pth.fn.bpod = 'Bpod.mat';
pth.fn.vid = 'Trajectories.mat';
pth.fn.obj =  ['data_structure_' anm '_' dt '.mat'];





















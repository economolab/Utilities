function showVid()


% parent = 'Z:\users\Mike\Projects\PTtagging\anm395710 - WR51\2017-10-26\vid';
% parent = 'Z:\users\Mike\Projects\PTtagging\anm395709 - WR52\2017-10-26\vid';
% parent = 'Z:\users\Mike\Projects\PTtagging\anm395712 - WR54\2017-10-26\vid';
% parent = 'Z:\users\Mike\Projects\PTtagging\anm395711 - WR53\2017-10-26\vid';

parent = 'C:\Users\labadmin\Desktop\DeepLabCut\Projects\Sydney\TrainingData';
stillDir = 'C:\Users\labadmin\Desktop\DeepLabCut\Projects\Sydney\TrainingData';

% parent = 'C:\Users\economom\Dropbox (HHMI)\BehaviorVideo\MNE\WR51';
% stillDir = 'C:\Users\economom\Dropbox (HHMI)\BehaviorVideo\MNE\StillImagesForTraining_below';

% parent = 'C:\Users\economom\Desktop';
% parent = 'Z:\users\tw\camera\TWdata\wr6\24-10-2017\0';

wsfn  = '*.h5';
% token = 'camera0';
token = 'cam_0_';

initFig(parent, wsfn, token, stillDir);



function initFig(parent, wsfn, token, stillDir)

h.parent = parent;
h.stillDir = stillDir;

[h.wsTimeStamp, h.wsTrialNum, h.bitcode, h.masking] = getTrialInfo(h.parent, wsfn);
[h.vidFns, h.vidTimeStamps] = findVideoFiles(h.parent, token);

if ~isempty(h.wsTrialNum)
    h.vidTrialNum = h.wsTrialNum(1:numel(h.vidFns));
else
    h.vidTrialNum = -1*ones(size(h.vidFns));
    h.wsTimeStamp = zeros(size(h.vidFns));
    h.wsTrialNum = -1*ones(size(h.vidFns));
end


bcol = [1 1 1];
fsize = 10;
h.frameRate = 400;

h.fig = figure;
set(h.fig, 'Units', 'Normalized', 'Position', [0.05 0.1 0.9 0.8], 'Color', bcol);
set(h.fig, 'Renderer', 'openGL');

h.ax(1) = axes;
set(h.ax(1), 'Units', 'Normalized', 'Position', [0 0.1 0.6 0.9]);
set(h.fig, 'WindowScrollWheelFcn', {@figScroll, h.fig});

h.vid = zeros(100, 100, 3);
h.medframe = mean(h.vid,3);
h.vidplot = imagesc(h.vid(:,:,1));
colormap(gray);
axis off;
axis image;


h.ax(2) = axes;
set(h.ax(2), 'Units', 'Normalized', 'Position', [0.6 0.6 0.28 0.38]);
set(h.ax(2), 'YTick', [], 'YTickLabel', [])

h.ax(3) = axes;
set(h.ax(3), 'Units', 'Normalized', 'Position', [0.6 0.15 0.28 0.24]);
set(h.ax(3), 'YTick', [], 'YTickLabel', [])


num = 1;
h.rect(num) = imrect(h.ax(num), [1 1 10 10]);
setColor(h.rect(num), 'g');
addNewPositionCallback(h.rect,@(x) moveROI(h.fig));


frames = size(h.vid,3);
h.frame_slider = uicontrol('Style', 'Slider', 'Units', 'normalized', 'Position', ...
    [0.13 0.05 0.67 0.03], 'Min', 1, 'Max', frames, 'Value', 1, 'SliderStep', 1./(max(frames-1, 1))+[0 0], ...
    'Callback', {@sliderCallback,gcf}, 'BackgroundColor', [0.75 0.75 0.75]);

uicontrol('Style', 'text', 'String', 'Frame: ','Units','normalized' ...
    ,'BackgroundColor',bcol,'Position',[0.01 0.05 0.05 0.03],'HorizontalAlignment','Right', ...
    'FontSize', fsize);


str = cell(numel(h.vidFns), 1);
for i = 1:numel(str)
   str{i} = [num2str(i) ': ' h.vidFns{i}];
end

h.framenumtext = uicontrol('Style', 'text', 'String', '1/3 ','Units','normalized' ...
    ,'BackgroundColor',bcol,'Position',[0.065 0.05 0.06 0.03],'HorizontalAlignment','Left', ...
    'FontSize', fsize);

h.file_listbox = uicontrol('Style', 'listbox', 'Units', 'normalized', 'Position', ...
    [0.9 0.06 0.1 0.9], 'String', str, 'Value', 1, 'BackgroundColor', [1 1 1], ...
    'Max', 3, 'Min', 1);

uicontrol('Style', 'text', 'String', 'Files','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.9 0.97 0.1 0.02],'HorizontalAlignment','Center', ...
    'FontSize', 10, 'FontWeight', 'Bold');

h.load_button = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.9 0.01 0.1 0.04], 'String', 'Load Vid','FontWeight','Bold', 'Callback', ...
    {@load_button_fcn,gcf});

h.medfilt_edit = uicontrol('Style', 'Edit', 'String', '3', 'Units', 'Normalized', ...
    'BackgroundColor', [1 1 1], 'Position', [0.65 0.52 0.05 0.025], 'HorizontalAlignment', 'Left', ...
    'Callback', {@updateTsAx,gcf});
uicontrol('Style', 'text', 'String', 'Med filt: ','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.6 0.52 0.05 0.025],'HorizontalAlignment','Right', ...
    'FontSize', fsize, 'FontWeight', 'Bold');

h.scrollSpeed_edit = uicontrol('Style', 'Edit', 'String', '10', 'Units', 'Normalized', ...
    'BackgroundColor', [1 1 1], 'Position', [0.65 0.49 0.05 0.025], 'HorizontalAlignment', 'Left');
uicontrol('Style', 'text', 'String', 'Scroll speed: ','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.6 0.49 0.05 0.025],'HorizontalAlignment','Right', ...
    'FontSize', fsize, 'FontWeight', 'Bold');

h.trialShift_edit = uicontrol('Style', 'Edit', 'String', '0', 'Units', 'Normalized', ...
    'BackgroundColor', [1 1 1], 'Position', [0.65 0.46 0.05 0.025], 'HorizontalAlignment', 'Left', ...
    'Callback', {@updateTrialPlot,gcf});
uicontrol('Style', 'text', 'String', 'Trial shift: ','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.6 0.46 0.05 0.025],'HorizontalAlignment','Right', ...
    'FontSize', fsize, 'FontWeight', 'Bold');

h.diffImage_checkbox = uicontrol('Style', 'checkbox',  'Units', 'normalized', 'Position', ...
    [0.6 0.42 0.15 0.025], 'String', 'Difference image', 'Backgroundcolor', bcol, 'FontSize', fsize, ...
    'Callback', {@updateVidAx,gcf});

h.saveROI_button = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.725 0.4 0.15 0.04], 'String', 'Save ROI','FontWeight','Bold', 'Callback', ...
    {@saveROI_button_fcn,gcf});

h.saveFrame_button = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.725 0.45 0.15 0.04], 'String', 'Save still image','FontWeight','Bold', 'Callback', ...
    {@saveFrame_button_fcn,gcf});

h.saveFrameInterval_button = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.725 0.5 0.15 0.04], 'String', 'Save by Interval','FontWeight','Bold', 'Callback', ...
    {@saveFrameInterval_button_fcn,gcf});

h.saveFrameROI_button = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.725 0.55 0.15 0.04], 'String', 'Save by ROI','FontWeight','Bold', 'Callback', ...
    {@saveFrameROI_button_fcn,gcf});


guidata(h.fig, h);
updateTrialPlot([], [], h.fig)





function saveFrameROI_button_fcn(~, ~, fig)

h = guidata(fig);


roiPos = getPosition(h.rect);
xix = round(roiPos(1):roiPos(1)+roiPos(3));
yix = round(roiPos(2):roiPos(2)+roiPos(4));

dat = bsxfun(@minus, h.vid(yix,xix,:), h.medframe(yix, xix));

ts = squeeze(mean(mean(abs(dat), 1), 2));
mf = round(str2double(get(h.medfilt_edit, 'String')));

if mf>0
    ts = medfilt1(ts, mf);
end

prompt = {'# Frames:'};
title = 'Input';
dims = [1 20];
definput = {'25'};
answer = inputdlg(prompt,title,dims,definput);
Npick = str2double(answer);


vals = linspace(min(ts), max(ts), Npick);
frames = zeros(Npick, 1);
for i = 1:Npick
    [~, frames(i)] = min(abs(ts-vals(i)));
end
frames = unique(frames);  

figure; plot(ts); hold on; plot(frames, ts(frames), 'r.');


for i = 1:numel(frames)
    
    vidval = get(h.file_listbox, 'value');
    vidfn = h.vidFns{vidval};
    
    fileparts = strsplit(h.parent, filesep);
    tokens = {'WR'; '201'; 'anm'};
    
    fn = [vidfn(1:end-4) '-frame' num2str(frames(i)) '.png'];
    for j = 1:numel(tokens)
        if sum(~cellfun(@isempty, strfind(fileparts, tokens{j})))>0
            str = fileparts{~cellfun(@isempty, strfind(fileparts, tokens{j}))};
            fn = [str '_' fn];
        end
    end
    
    % dirname = fullfile(h.parent, 'stills');
    dirname = h.stillDir;
    if ~isdir(dirname)
        mkdir(dirname);
    end
    
    if get(h.diffImage_checkbox, 'Value')
        im = abs(h.vid(:,:,frames(i)) - h.medframe);
        im = uint8(floor(255.*(im - min(im(:)))./(max(im(:)) - min(im(:)))));
    else
        im = uint8(h.vid(:,:,frames(i)));
    end
    rgb = im;
    rgb(:, :, 2) = im;
    rgb(:, :, 3) = im;
    % figure(324); imagesc(im); colorbar; colormap(gray);
    imwrite(rgb,fullfile(dirname, fn));
end






function saveFrameInterval_button_fcn(~, ~, fig)

h = guidata(fig);

curframe = round(get(h.frame_slider, 'value'));

prompt = {'Frames start:','Frames end:', 'Increment:'};
title = 'Input';
dims = [1 20];
definput = {num2str(curframe),num2str(curframe+100),'2'};
answer = inputdlg(prompt,title,dims,definput);
frames = str2double(answer{1}):str2double(answer{3}):str2double(answer{2});



for i = 1:numel(frames)
    
    vidval = get(h.file_listbox, 'value');
    vidfn = h.vidFns{vidval};
    
    fileparts = strsplit(h.parent, filesep);
    tokens = {'WR'; '201'; 'anm'};
    
    fn = [vidfn(1:end-4) '-frame' num2str(frames(i)) '.png'];
    for j = 1:numel(tokens)
        if sum(~cellfun(@isempty, strfind(fileparts, tokens{j})))>0
            str = fileparts{~cellfun(@isempty, strfind(fileparts, tokens{j}))};
            fn = [str '_' fn];
        end
    end
    
    % dirname = fullfile(h.parent, 'stills');
    dirname = h.stillDir;
    if ~isdir(dirname)
        mkdir(dirname);
    end
    
    if get(h.diffImage_checkbox, 'Value')
        im = abs(h.vid(:,:,frames(i)) - h.medframe);
        im = uint8(floor(255.*(im - min(im(:)))./(max(im(:)) - min(im(:)))));
    else
        im = uint8(h.vid(:,:,frames(i)));
    end
    rgb = im;
    rgb(:, :, 2) = im;
    rgb(:, :, 3) = im;
    % figure(324); imagesc(im); colorbar; colormap(gray);
    imwrite(rgb,fullfile(dirname, fn));
end









function saveFrame_button_fcn(~, ~, fig)

h = guidata(fig);

vidval = get(h.file_listbox, 'value');
vidfn = h.vidFns{vidval};

frame = round(get(h.frame_slider, 'value'));
fileparts = strsplit(h.parent, filesep);
tokens = {'WR'; '201'; 'anm'};

fn = [vidfn(1:end-4) '-frame' num2str(frame) '.png'];
for i = 1:numel(tokens)
    if sum(~cellfun(@isempty, strfind(fileparts, tokens{i})))>0
        str = fileparts{~cellfun(@isempty, strfind(fileparts, tokens{i}))};
        fn = [str '_' fn];
    end
end

% dirname = fullfile(h.parent, 'stills');
dirname = h.stillDir;
if ~isdir(dirname)
    mkdir(dirname);
end

if get(h.diffImage_checkbox, 'Value')
    im = abs(h.vid(:,:,frame) - h.medframe);
else
    im = h.vid(:,:,frame);
end
im = uint8(floor(255.*(im - min(im(:)))./(max(im(:)) - min(im(:)))));
rgb = im;
rgb(:, :, 2) = im;
rgb(:, :, 3) = im;
% figure(324); imagesc(im); colorbar; colormap(gray);
imwrite(rgb,fullfile(dirname, fn));







function saveROI_button_fcn(~, ~, fig)

h = guidata(fig);
% 
% v = VideoWriter(fullfile(h.parent, 'newfile_side_stim2.avi'));
% v.FrameRate = 333;
% 
% open(v);
% dat = uint8(h.vid);
% dat = permute(dat, [1 2 4 3]);
% writeVideo(v,dat)
% close(v)


roiPos = getPosition(h.rect);
ix.x = round(roiPos(1):roiPos(1)+roiPos(3));
ix.y = round(roiPos(2):roiPos(2)+roiPos(4));
meanImage = h.medframe;
parent = h.parent;

wsTimeStamp = h.wsTimeStamp;
wsTrialNum  = h.wsTrialNum;
vidTimeStamps = h.vidTimeStamps;
vidFns = h.vidFns;
wsBitcode = h.bitcode;
vidTrialNum = h.vidTrialNum;

fn = fullfile(h.parent, 'ROI.mat');
save(fn, 'wsTimeStamp', 'wsTrialNum', 'wsBitcode', 'vidTimeStamps', 'vidFns', 'vidTrialNum', 'parent', 'ix', 'meanImage');
disp(['Wrote ' fn]);





function moveROI(fig)

updateTsAx([], [], fig);




function sliderCallback(~, ~, fig)
updateVidAx([], [], fig);
updateTsAx([], [], fig);






function updateTrialPlot(~, ~, fig)
h = guidata(fig);

shift = str2double(get(h.trialShift_edit, 'String'));
% hold(h.ax(3), 'off');
% plot(h.ax(3), 1:numel(h.wsTimeStamp)-1, diff(h.wsTimeStamp));
% 
% hold(h.ax(3), 'on');
% plot(h.ax(3), shift+(1:numel(h.vidTimeStamps)-1), diff(h.vidTimeStamps));


hold(h.ax(3), 'off');
N = numel(h.vidTimeStamps);

if shift+N>numel(h.wsTimeStamp)
    disp('Invalid shift');
    set(h.trialShift_edit, 'String', '0');
    shift = 0;
end


plot(h.ax(3), h.vidTimeStamps - h.wsTimeStamp(shift+(1:N)));
axis(h.ax(3), 'tight');

str = cell(numel(h.vidFns), 1);
for i = 1:N
   str{i} = [num2str(h.wsTrialNum(shift+i)) ': ' h.vidFns{i}];
end
h.vidTrialNum = h.wsTrialNum(shift+1:N);
set(h.file_listbox, 'String', str);

guidata(fig, h);





function updateTsAx(~, ~, fig)
h = guidata(fig);


roiPos = getPosition(h.rect);
xix = round(roiPos(1):roiPos(1)+roiPos(3));
yix = round(roiPos(2):roiPos(2)+roiPos(4));

dat = bsxfun(@minus, h.vid(yix,xix,:), h.medframe(yix, xix));

ts = squeeze(mean(mean(abs(dat), 1), 2));
mf = round(str2double(get(h.medfilt_edit, 'String')));

if mf>0
    ts = medfilt1(ts, mf);
    frame = round(get(h.frame_slider, 'value'));
end

hold(h.ax(2), 'off');
plot(h.ax(2), (1:numel(ts))./h.frameRate, ts, 'k-');
yl = [min(ts) max(ts)];
hold(h.ax(2), 'on');
plot(h.ax(2), [frame frame]./h.frameRate, yl, 'r-');






function figScroll(~,evnt,fig)

h = guidata(fig);


frame = round(get(h.frame_slider, 'value'));
frame = frame+str2double(get(h.scrollSpeed_edit, 'String')).*evnt.VerticalScrollCount;
    maxframe = size(h.vid,3);

if frame<1
    frame = 1;
end
if frame>maxframe
    frame = maxframe;
end

set(h.frame_slider, 'value', frame);

updateVidAx([], [], fig);
updateTsAx([], [], fig);





function load_button_fcn(~, ~, fig)

h = guidata(fig);

vidstr = get(h.file_listbox, 'string');
vidval = get(h.file_listbox, 'value');

vidfn = fullfile(h.parent, h.vidFns{vidval});

h.vidObj = VideoReader(vidfn);
frames = h.vidObj.Duration*h.vidObj.FrameRate;

h.vid = zeros(h.vidObj.Height, h.vidObj.Width, frames);

for i = 1:frames
    h.vid(:,:,i) = mean(readFrame(h.vidObj), 3);
end

set(h.frame_slider, 'Max', size(h.vid, 3), 'SliderStep', [1 1]./size(h.vid,3));
h.medframe = median(h.vid, 3);

guidata(fig, h);
updateTsAx([], [], fig);
updateVidAx([], [], fig);




function updateVidAx(~, ~, fig)
h = guidata(fig);

frame = round(get(h.frame_slider, 'value'));

if get(h.diffImage_checkbox, 'Value')
    set(h.vidplot, 'Cdata', abs(h.vid(:,:,frame) - h.medframe));
    set(h.ax(1), 'Clim', [-100 100]);
    colormap(hsv);
else
    set(h.vidplot, 'Cdata', h.vid(:,:,frame));
    set(h.ax(1), 'Clim', [0 255]);
    colormap(gray);
end
set(h.framenumtext, 'String', [num2str(frame) '/' num2str(size(h.vid, 3))]);





function [vidFns, vidTimeStamps] = findVideoFiles(parent, token)

contents = dir(parent);
good = false(numel(contents), 1);
for i = 1:numel(good)
   if ~isempty(strfind(contents(i).name, token));
       good(i) = true;
   end
end


contents = contents(good);
vidTimeStamps = zeros(numel(contents), 1);
vidFns = cell(numel(contents), 1);
for i = 1:numel(contents)
    vidFns{i} = contents(i).name;
    vidTimeStamps(i) = contents(i).datenum;
end

vidTimeStamps = vidTimeStamps*60*60*24;





function [timeStamp, trialNum, bitcode, masking] = getTrialInfo(parent, wsfn)

datfile = ls(fullfile(parent, wsfn));

if size(datfile, 1)~=1
    timeStamp = [];
    trialNum = [];
    bitcode = [];
    masking = [];
    return;
end
dat=ws.loadDataFile(fullfile(parent, datfile));

fnames = fieldnames(dat);

dataFields = ~cellfun(@isempty, (strfind(fnames, 'sweep')));
N = sum(dataFields);

timeStamp = -1.*ones(N, 1);
trialNum = -1.*ones(N, 1);
bitcode = cell(sum(dataFields), 1);
masking = cell(sum(dataFields), 1);
for i = find(dataFields')
    
    timeStamp(i) = dat.(fnames{i}).timestamp;
    bitcode{i} = dat.(fnames{i}).analogScans(:,1);
    masking{i} = dat.(fnames{i}).analogScans(:,4);
    
    [trialNum(i), ~] = getBitCode(bitcode{i}, dat.header.Acquisition.SampleRate);
end
timeStamp = timeStamp + 24*60*60*datenum(dat.header.ClockAtRunStart);






function [num, ix] = getBitCode(bitcode, fs)

try
    sampPerBit = round(fs*0.00715); %7 ms per bit
    Nbits = 10;
    
    threshold = 1;
    
    above = bitcode>threshold;
    above = [false; above; false];
    edges = diff(above);
    rising = find(edges==1);
    falling = find(edges==-1);
    width = falling-rising;
    landmark = find(width>1000, 1, 'first');
    ix = rising(landmark);
    
    offset = round(0.07224*fs);
    ix = ix-offset;
    a = bitcode(ix:ix-1+sampPerBit*Nbits);
    a = a-median(a);
    bitmat = reshape(a, sampPerBit, Nbits);

    bits = mean(bitmat, 1)>threshold;
    num = bin2dec(num2str(fliplr(bits)')');
    
%     figure(1); hold off; plot(bitcode(ix:ix+2000));
    disp(['Assigning solo trial #' num2str(num)]);
    
    ix = ix+offset;
        
catch
    disp('!!Could not determine trial number!!!  Assigning -1');
    num = -1;
    ix = 0;
end




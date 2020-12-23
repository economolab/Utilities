function loadVideoTrajectories(pth)

for i = 1:numel(pth.vid)
    fn = findFiles(pth.vid{i});
    traj{i} = loadCSVFeatTrajectories(pth.vid{i}, fn);
end

save(fullfile(pth.sv, pth.fn.vid), 'traj');


view = 1;
figure; hold on;
a = NaN.*zeros(5000, numel(traj{view}));
Nf = numel(traj{view});
for i = 1:numel(traj{view})
    feat = 2;
%     plot(traj{view}(i).ts(:, 1, feat), 'r'); 
clr = hsv2rgb([i/1000 1 1]);
    plot(traj{view}(i).ts(:, 2, feat), 'Color', clr); 
    a(1:size(traj{view}(i).ts, 1), i) = traj{view}(i).ts(:, 2, feat);
end
figure;
imagesc(a);





function traj = loadCSVFeatTrajectories(parent, fn)
pthresh = 0.9;


for i = 1:numel(fn)
    disp(['Loading ' fn{i}]);
    disp(['     File ' num2str(i) '/' num2str(numel(fn))])
    fid = fopen(fullfile(parent,fn{i}));
    h1 = fgetl(fid);
    h2 = fgetl(fid);
    h3 = fgetl(fid);
    cnt = 0;
    while(~feof(fid))
        str = fgetl(fid);
        cnt = cnt+1;
    end
    fclose(fid);
    
    newStr = split(h2,',');
    bodyparts = newStr(2:3:end);
    
    Nframes = cnt;
    Nparts = numel(bodyparts);
    ts = zeros(Nframes, 3, Nparts);
    
    fid = fopen(fullfile(parent,fn{i}));
    h1 = fgetl(fid);
    h2 = fgetl(fid);
    h3 = fgetl(fid);
    cnt = 0;
    while(~feof(fid))
        str = fgetl(fid);
        cnt = cnt+1;
        
        dat = cellfun(@str2double, split(str, ','));
        ts(cnt, :, :) = reshape(dat(2:end)', 3, Nparts);
        
    end
    fclose(fid);
    
    for j = 1:Nparts
        nans = ts(:, 3, j)<pthresh;
        ts(nans, 1, j) = NaN;
        ts(nans, 2, j) = NaN;
    end
    
    traj(i).fn = fn{i};
    traj(i).featNames = bodyparts;
    traj(i).ts = ts;
    
%     if i==numel(fn)
%         'hi'
%     end
end



% 
% figure; hold on;
% a = NaN.*zeros(5000, numel(traj));
% for i = 1:numel(traj)
%     feat = 8;
%     plot(traj(i).ts(:, 1, feat), 'r'); 
%     plot(traj(i).ts(:, 2, feat), 'b'); 
%     a(1:size(traj(i).ts, 1), i) = traj(i).ts(:, 2, feat);
% end
% figure;
% imagesc(a);



function traj = loadH5FeatTrajectories(parent, fn)
NptsPerFeat = 3;


for i = 1:numel(fn)
    disp(['Loading ' fn{i}]);
    info = hdf5info(fullfile(parent,fn{i}));
    data = hdf5read(info.GroupHierarchy.Groups(1).Datasets);
    
    Nframes = numel(data);
    traj(i).ts = cell(Nfeat, 1);
    for j = 1:Nfeat
        traj(i).fn = fn{i};
        traj(i).featNames = featNames;
        traj(i).ts{j} = zeros(Nframes, NptsPerFeat);
    end
    
    for j = 1:Nframes
        for k = 1:Nfeat
            try
                traj(i).ts{k}(j, :) = (data(j).Data{2}.Data(1+(k-1)*NptsPerFeat:k*NptsPerFeat))';
            catch
                'hi'
            end
        end
    end
end





function fn = findFiles(parent)


contents = dir(fullfile(parent, '*.csv'));

fn = cell(numel(contents), 1);
for i = 1:numel(contents)
    fn{i} = contents(i).name;
end










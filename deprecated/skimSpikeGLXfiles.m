function skimSpikeGLXfiles()
pth = 'Y:\EEL\Experiments\EEL2\SpikeGLX\2019-09-10';


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

events = zeros(numel(contents), 1);
for i = 1:numel(contents)-1
    events(i+1) = getMetaData(fullfile(pth, contents(i).name), 'fileTimeSecs');
end
events = cumsum(events);
save(fullfile(pth, 'events.mat'), 'events');



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




function out = getDataStructures(pth)
% input path to Analysis folder
% Analysis folder should contain date folders in format YEAR-MO-DA
data_dirs = dir(pth);
for i=1:numel(data_dirs)
   if length(data_dirs(i).name) == 10
       data_files = dir(fullfile(pth, data_dirs(i).name));
       for j = 1:numel(data_files)
           if contains(data_files(j).name, 'data_structure')
               
               fn{i} = fullfile(pth, data_dirs(i).name, data_files(j).name);
               date{i} = data_dirs(i).name;
           end
       end
   end
end
fn = fn(~cellfun(@isempty, fn));
date = date(~cellfun(@isempty, date));

for i=1:numel(fn)
   a(i) = load(fn{i}); 
   a(i).obj.date = date{i};
end

for i = 1:numel(a)
    out(i) = a(i).obj;
end
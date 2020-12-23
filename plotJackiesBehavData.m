parent = 'Z:\JEB\Experiments\JEB2\Analysis';
animal = 'JEB2';
dates = {'2020-10-09','2020-10-10','2020-10-11','2020-10-12','2020-10-17','2020-10-18'};

contents = dir(parent);

%Load all file names from the Analysis Folder
fn = {};
for i = 1:numel(contents)
    if exist(fullfile(parent, contents(i).name), 'dir') && numel(contents(i).name)>2
        fn{end+1} = fullfile(parent, contents(i).name, 'bpod.mat');
    end
end
%Change depending on which dates you want to analyze
fn = fn(end-5:end);

%Extract the bp data from each file 
for i = 1:numel(fn)
    a = load(fn{i});
    bp(i) = a.bp;
end

Nshow = 10;  %How many trials after a stim you want to look at
Nsess = numel(bp);
pct.R = zeros(Nsess,Nshow+1);
pct.L = zeros(Nsess,Nshow+1);

for i = 1:Nsess
    
    b = bp(i);
    
    b.stim.sinceStim = getLastStim(b);
    
    for j = 0:Nshow
        f = b.stim.sinceStim == j;
        pct.R(i,j+1) = sum(~b.no&b.R&b.hit&f)/sum(~b.no&b.R&f);
        pct.L(i,j+1) = sum(~b.no&b.L&b.hit&f)/sum(~b.no&b.L&f);
    end
    
    
%     ix = b.stim.sinceStim>Nshow;
%     
%     pct.R = sum(b.no&b.R&ix)./sum(b.R&ix);
%     pct.L = sum(b.no&b.L&ix)./sum(b.L&ix);
% 
%     figure(1);
%     plot(0, pct.L, 'ro'); hold on;  plot(pct.R, 'bo');
%     
%     for j = 1:Nshow
%         ix = (b.stim.sinceStim+1)==j;
%         stimpct(j, i, 1) = sum(b.no&b.R&ix)./sum(b.R&ix);
%         stimpct(j, i, 2) = sum(b.no&b.L&ix)./sum(b.L&ix);
%     end
end
% 
figure(2);
for w = 1:Nsess
    plot(0:Nshow,pct.R(w,:),'r:'); hold on; plot(0:Nshow,pct.L(w,:),'b:'); hold on;
end
legend('Right','Left')
plot(0:Nshow,mean(pct.R(2:4,:),1),'r'); hold on; plot(0:Nshow,mean(pct.L(2:4,:),1),'b');
hold off;
xlabel('# trials since stim');
ylabel('Behav Performance')
title([animal,' ', dates])

% plot(nanmean(stimpct(:,:,1), 2), 'b.-');
% plot(nanmean(stimpct(:,:,2), 2), 'r.-');
% 
% xlim([-1 Nshow+1]);

function sinceStim = getLastStim(b)

%Find the indices of all stim trials
stimnums = find(b.stim.enable);

sinceStim = -inf*ones(b.Ntrials, 1);

%For each trial...
for j = 1:b.Ntrials
    %If the trial occurred later than the first stim trial
    if j>min(stimnums)
        %Get the distance between trial and all of the stim trial indices
        tmp = j-stimnums;
        %Discard any negative distances
        tmp(tmp<0) = inf;
        %Will give you a distance since the last stim trial
        sinceStim(j) = min(tmp);
    end
end
end
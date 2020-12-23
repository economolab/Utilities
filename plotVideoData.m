function plotVideoData()

mouse = 'EEL6';
pth = ['E:\', mouse, '\Analysis'];

obj = getDataStructures(pth);

'hi'
for i = 1:3%numel(obj)
    thisobj = obj(i);
    R = thisobj.bp.R;
    L = thisobj.bp.L;
    hit = thisobj.bp.hit;
    early = thisobj.bp.early;
    lickL = thisobj.bp.ev.lickL;
    lickR = thisobj.bp.ev.lickR;
    stim = thisobj.bp.stim;
    no = thisobj.bp.no;
    
    med = 7;
    view =2;
    % side: view = 1
    % bottom: view = 2
    
    bodyparts = 1:3;
    offset = 25;
    figure;
    plotStackedTrajectories(thisobj, bodyparts(1), view, {L&stim.num==1; L&stim.num==2; L&stim.num==3; L&stim.num==4}, med, {[1 0 0]; [0 0 1]; [0 0 0]; [0 0.75 0]}, offset);
end
'hi'


figure;
plotTrajectories(thisobj, bodyparts, view, R&hit, med, 'b', 1);
plotTrajectories(thisobj, bodyparts, view, L&hit, med, 'r', 1);
'hi'

% %view 2 bodyparts
%     1{'top_tongue'       }
%     2{'topleft_tongue'   }
%     3{'bottom_tongue'    }
%     4{'leftbottom_tongue'}
%     5{'top_paw'          }
%     6{'bottom_paw'       }
%     7{'lickport'         }
%     8{'jaw'              }

% single day
figure;
bodyparts = [1];
view = 2;
plotTrajectories(thisobj, bodyparts, view, L&hit&stim.num==0, med, 'r', 0.75);
plotTrajectories(thisobj, bodyparts, view, L&stim.num==1, med, 'g', 1.5);
plotTrajectories(thisobj, bodyparts, view, L&stim.num==2, med, 'm', 1.5);
plotTrajectories(thisobj, bodyparts, view, L&stim.num==3, med, 'c', 1.5);
plotTrajectories(thisobj, bodyparts, view, L&stim.num==4, med, 'y', 1.5);
figTitle = [thisobj.traj{view}(1).featNames{bodyparts(1)}];
for i = 2:numel(bodyparts)
    figTitle = [figTitle, ', ', thisobj.traj{view}(1).featNames{bodyparts(i+1)}];
end
title(figTitle)


% multiday
ix = 1:5;
view = 2;
med = 7;
bodyparts = 1;
plotStackedMultidayTraj(obj, ix, bodyparts, view, 'L', med, 'r', 1.5);
bodyparts = 3;
plotStackedMultidayTraj(obj, ix, bodyparts, view, 'R', med, 'b', 1.5);


% plot by stim type
ix = 1:5;
view = 2;
med = 7;
bodyparts = 1;
plotStimtypesSubplot(obj, ix, bodyparts, view, 'L', med, 'r', 1.5)
sgtitle([mouse, ' Left Trials'], 'FontWeight', 'bold');
bodyparts = 3;
plotStimtypesSubplot(obj, ix, bodyparts, view, 'R', med, 'b', 1.5)
sgtitle([mouse, ' Right Trials'], 'FontWeight', 'bold');


% bilateral vs unilateral
med = 7;
view = 2;
bodyparts = 1;
ix1 = 1:3;
ix2 = 4:5;
stimNum = '3';
plotHemisphereMultiday(obj, ix1, ix2, bodyparts, view, 'L', stimNum, med, 'r', 1.5)
sgtitle([mouse, ' Stim ', stimNum], 'FontWeight', 'bold');
bodyparts = 3;
plotHemisphereMultiday(obj, ix1, ix2, bodyparts, view, 'R', stimNum, med, 'b', 1.5)
sgtitle([mouse, ' Stim ', stimNum], 'FontWeight', 'bold');


figure;
bodyparts = 1;
plotTrajMultiday(obj, 1:5, bodyparts, view, 'L&hit&stim.num==0', med, 'r', 0.5);
plotTrajMultiday(obj, 1:3, bodyparts, view, 'L&stim.num==3', med, [1 0 1], 2.5);
plotTrajMultiday(obj, 4:5, bodyparts, view, 'L&stim.num==3', med, [1 1 0], 2.5);

bodyparts = [3 5 6];
for i=1:4
    figure;
    plotTrajMultiday(obj, 1:5, i, view, 'R&hit&stim.num==0', med, 'b', 0.5);
    plotTrajMultiday(obj, 1:3, i, view, 'R&hit&stim.num==3', med, [1 0 1], 2.5);
    plotTrajMultiday(obj, 4:5, i, view, 'R&hit&stim.num==3', med, [0 1 1], 2.5);
    
    figure;
    plotTrajMultiday(obj, 1:5, i, view, 'L&hit&stim.num==0', med, 'r', 0.5);
    plotTrajMultiday(obj, 1:3, i, view, 'L&stim.num==3', med, [1 0 1], 2.5);
    plotTrajMultiday(obj, 4:5, i, view, 'L&stim.num==3', med, [1 1 0], 2.5);
end
'hi'












function plotHemisphereMultiday(obj, ix1, ix2, bodyparts, view, whichTT, stimNum, med, clr, lineWid)
figure;
subplot(1,2,1)
plotTrajMultiday(obj, ix1, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix1, bodyparts, view, [whichTT, '&hit&stim.num==', stimNum], med, [1 0 1], lineWid);
title('Bilateral');
subplot(1,2,2)
plotTrajMultiday(obj, ix2, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix2, bodyparts, view, [whichTT, '&hit&stim.num==', stimNum], med, [0 1 1], lineWid);
title('Unilateral');



function plotStackedMultidayTraj(obj, ix, bodyparts, view, whichTT, med, clr, lineWid)
figure;
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==1'], med, 'c', lineWid);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==2'], med, 'm', lineWid);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==3'], med, 'g', lineWid);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==4'], med, 'y', lineWid);


function plotStimtypesSubplot(obj, ix, bodyparts, view, whichTT, med, clr, lineWid)

figure;
subplot(2,2,1);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&stim.num==1'], med, 'c', lineWid);
title('Stim 1')
subplot(2,2,2);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&stim.num==2'], med, 'm', lineWid);
title('Stim 2')
subplot(2,2,3);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==3'], med, 'g', lineWid);
title('Stim 3')
% plotTrajMultiday(obj, ix, bodyparts, view, ['L&~hit&stim.num==3'], med, 'b', 1.5);
subplot(2,2,4);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&hit&stim.num==0'], med, clr, 0.5);
plotTrajMultiday(obj, ix, bodyparts, view, [whichTT, '&stim.num==4'], med, 'y', lineWid);
title('Stim 4');
sgtitle([whichTT, ' Trials']);




function plotTrajMultiday(obj, ix, bodyparts, view, whichTT, med, clr, lineWid)
hold on;

for i=ix
    thisobj = obj(i);
    R = thisobj.bp.R;
    L = thisobj.bp.L;
    hit = thisobj.bp.hit;
    early = thisobj.bp.early;
    lickL = thisobj.bp.ev.lickL;
    lickR = thisobj.bp.ev.lickR;
    stim = thisobj.bp.stim;
    no = thisobj.bp.no;
    
    thisTT = eval(whichTT);
    
    plotTrajectories(thisobj, bodyparts, view, thisTT, med, clr, lineWid);
end
    



function plotStackedTrajectories(obj, bodypart, view, trials, med, clrs, offset)

Ntrials = numel(obj.traj{1});

dt = 1./400;
yoffset = 0;
for k = 1:2
    subplot(1,2,k); hold on;
    for j = 1:numel(trials)
        for i = 1:Ntrials
            if ~any(trials{j}(i))
                continue;
            end
            
            ts = medfilt1(obj.traj{view}(i).ts(:, k, bodypart), med, [], 1);
            time = dt.*(1:numel(ts));
            plot(time, ts+yoffset*offset, 'Color', clrs{j});
            %     xoffset = xoffset+time(end);
            yoffset = yoffset+1;
        end
    end
end
% xlim([0 5]);
yl = ylim();
plot(2.48+[0 0], yl, 'k:');



function ts = getTrajectories(obj, bodypart, view, trials)

Ntrials = sum(trials);
len = 0;
for i = 1:numel(obj.traj{view})
    len = max(len, size(obj.traj{view}(i).ts, 1));
end

ts = NaN.*zeros(len, 3, Ntrials);
cnt = 0;

trialnum = find(trials);
for i = 1:Ntrials
    
    cnt = cnt+1;
    dat = obj.traj{view}(trialnum(i)).ts(:, :, bodypart);
    ts(1:size(dat, 1), :, cnt) = dat;
end


function ts = plotTrajectories(obj, bodyparts, view, trials, med, clr, lineWid)


for i = bodyparts
    bodypart = i;
    
    hold on;
    
    ts = medfilt1(getTrajectories(obj, bodypart, view, trials), med, [], 1);
    plot(squeeze(ts(:, 1, :)), -squeeze(ts(:, 2, :)), 'Color', clr, ...
        'LineWidth', lineWid)

end
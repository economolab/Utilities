% for transformatting
% deeplabcut labels to APT


clear ndata


% parent = ('C:\Users\labadmin\Desktop\DeepLabCut\Projects\Sangwook_v1\Sangwook_v1-SP-2020-02-12\videos\iteration5\EEL6_2020-02-18__cam_1_date_2020_02_18_time_14_09_26_v001DLC_resnet50_Sangwook_v1Feb12shuffle1_1030000.csv');
 parent = ('C:\Users\labadmin\Desktop\DeepLabCut\Projects\Sydney\Sydney-SMB-2019-12-04\videos\EEL8_2019-11-22__cam_0_date_2019_11_22_time_15_04_44_v001DLC_resnet50_SydneyDec4shuffle1_1030000.csv');
[ndata, text, alldata] = xlsread(parent);

if size(alldata,2) == 25
for i = 1:size(ndata,1)
    for n = 4:3:25
    if ndata(i,n)<0.95
        ndata(i,n-1)=NaN;
        ndata(i,n-2)=NaN;
    end
    
    end
end

for k = 1:8
    ndata(:,28-3*k)=[];
end

ndata(:,1)=[];
elseif size(alldata,2) == 13
    for i = 1:size(ndata,1)
    for n = 4:3:13
    if ndata(i,n)<0.95
        ndata(i,n-1)=NaN;
        ndata(i,n-2)=NaN;
    end
    
    end
    end

for k = 1:4
    ndata(:,16-3*k)=[];
end

ndata(:,1)=[];
else 
    disp('this is not convertable because of different number of landmarks');
end


pTrkFrm = 1:size(ndata,1);% # frames
numframes = size(pTrkFrm,2);
pTrkiPt = round(1:(size(ndata,2)/2));% # landmarks
numlandmarks = size(pTrkiPt,2);
pTrkiTgt = 1; % # target
numtargets = size(pTrkiTgt,2);

e = 2*pTrkiTgt;
f = numframes*numtargets; 

pTrkTS = zeros(numlandmarks,numframes);
for a = size(ndata,1):-1:0
    c = (1:(size(pTrkiPt,2)))';
    c(:,1) = 1;
    d = 1:size(ndata,1);
    d(1,:) = now-a/(24*3600);
    pTrkTS(:,:) = c*d; % timestamps  
end


pTrk = zeros(numlandmarks, e,f);
pTrk(:, 1, :) = ndata(:, 1:2:end)';
pTrk(:, 2, :) = ndata(:, 2:2:end)';

pTrkTag = zeros(numlandmarks,numframes);
u = pTrk(:,1,:);
pTrkTag(:,:) = isnan(u);
pTrkTag = logical(pTrkTag); % # occulusion

newName = num2str(parent);
save('newName11.trk','pTrk','pTrkFrm','pTrkTS','pTrkiPt','pTrkiTgt','pTrkTag')



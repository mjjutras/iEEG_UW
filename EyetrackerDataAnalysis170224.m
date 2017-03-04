Fs = 500;
dt = 1/Fs;

subjdir = 'R:\obregon\NovelvsRepeatTask\subjects';
subjid = '5d37ff';
rec_day = 'd4';
cd(fullfile(subjdir,subjid,'data',rec_day));

d = dir('* Samples.txt');
format = [repmat('%s',1,9) '%*[^\n]'];

blk = 1; % go from block 1 to 5 (6 & 7 are test blocks)
fid = fopen(d(blk).name);
C = textscan(fid,format,'HeaderLines',41,'Delimiter','\t');

% column 1: Time (in microseconds)
% column 4: Trial messages (BLK1, S02I23.bmp, etc.)
% column 8: x position
% column 9: y position

% indx1 
indx1 = find(~cellfun(@isempty,strfind(C{4},'.bmp')));

isi_inds = find(~cellfun(@isempty,strfind(C{4},'ISI')),length(indx1)-1,'last');
indx2 = [isi_inds; find(~cellfun(@isempty,strfind(C{4},['END_BLK' num2str(blk)])))];

img = 1; % image 21 is repeat of image 1, etc.
samples = indx1(img)+1:indx2(img)-1;
xpos = str2double(C{8}(samples));
ypos = str2double(C{9}(samples));
times = 10^-6 * str2double(C{1}(samples)); % converted to seconds

%% plays a movie of the eyetracker data (not really essential just interesting)
xh = [0 2000]; yh = [0 1200];

figure

clear M1
for frlop = 1:9000
    clf
    scatter(str2double(C{8}((frlop*10-9:frlop*10)+1)),str2double(C{9}((frlop*10-9:frlop*10)+1)),'.b')
    xlim(xh);ylim(yh);
    M1(frlop)=getframe;
end

% movie2avi('test.avi','FPS',20)
% disp('Created file')

%%
time = str2double(C{1}(2:end));
xpos = str2double(C{8}(2:end));
ypos = str2double(C{9}(2:end));

time = time(~isnan(ypos));
xpos = xpos(~isnan(ypos));
ypos = ypos(~isnan(ypos));

figure
ax(1) = subplot(2,1,1);
plot(time(1:9999),diff(time(1:10000)))
ax(2) = subplot(2,1,2);
plot(time(1:10000),[xpos(1:10000) ypos(1:10000)],'o-')
linkaxes(ax,'x')

% time values are off, delays between successive samples aren't consistent,
% fix this
% true sampling rate is 500 Hz (or close enough that it works as an
% approximation) so delay between samples should be 2 ms

figure;hold on
scatter(time,time,'.b')
scatter(time(1):2000:time(end),time(1):2000:time(end),'or')

newtime = time(1):2000:time(end);
newx = nan(size(newtime));
newy = nan(size(newtime));
parfor k = 1:length(newtime)
    if abs(newtime(k)-time(ft_nearest(time,newtime(k))))<1000
        newx(k) = xpos(ft_nearest(time,newtime(k)));
        newy(k) = ypos(ft_nearest(time,newtime(k)));
    end
end

newx = inpaint_nans(newx,2);
newy = inpaint_nans(newy,2);

%%
fltord = 60;
lowpasfrq = 25; % tried 18, might be too low;  seems to perform well
nyqfrq = 500 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);

% mark segments where eyedata (y data here) drops to 0; also mark 20
% samples (40 ms) before and 25 samples (50 ms) after each marked segment
% with 2 to account for saccades
% note that after filtering and calculating velocity, the time marked here
% only accounts for part of the saccade most of the time. I'm leaving a bit
% behind to go into ClusterFix. we can remove it later if we want
y_bound = double(newy==0);
for k = 1:length(y_bound)
    if k~=1 && y_bound(k)==1 && y_bound(k-1)==0 % start of marked segment
        if k>=21
            y_bound(k-1:-1:k-20) = 2;
        else
            y_bound(1:k-1) = 2;
        end
    elseif y_bound(k)==1 && y_bound(k+1)==0 % end of marked segment
        if k<=length(y_bound)-26
            y_bound(k+1:k+25) = 2;
        else
            y_bound(k+1:length(y_bound)) = 2;
            break
        end
    end
end

newx = [newx(100:-1:1) newx newx(end:-1:end-99)]; %add buffer for filtering
newy = [newy(100:-1:1) newy newy(end:-1:end-99)];   %add buffer for filtering
xss = filtfilt(flt,1,newx);
yss = filtfilt(flt,1,newy);
xss = xss(101:end-100); %remove buffer for filtering
yss = yss(101:end-100); %remove buffer for filtering
newx = newx(101:end-100);
newy = newy(101:end-100);

figure
ax(1) = subplot(2,1,1);
plot(newtime(1:10000),[newx(1:10000); newy(1:10000)],'o-')
ax(2) = subplot(2,1,2);
plot(newtime(1:10000),[xss(1:10000); yss(1:10000)],'o-')
linkaxes(ax,'x')

% velx = diff(xss);
% vely = diff(yss);
% vel = sqrt(velx.^2+vely.^2);

% transform all the variables to remove the "boundary" data
xss(logical(y_bound)) = nan;
yss(logical(y_bound)) = nan;
velx = diff(xss);
vely = diff(yss);
vel = sqrt(velx.^2+vely.^2);

figure
ax(1) = subplot(2,1,1);
plot(newtime(1:30000),[xss(1:30000); yss(1:30000)])
ax(2) = subplot(2,1,2);
plot(newtime(1:30000),vel(1:30000))
linkaxes(ax,'x')

accel = abs(diff(vel));
angle = 180*atan2(vely,velx)/pi;
%     vel = vel(200:end-200); %we do not care about the 1st or last 200 ms, line 208 assumes this too!
vel=vel(1:end-1);
%     accel = accel(200:end-199);
rot = zeros(1,length(xss)-2);
dist = zeros(1,length(xss)-2);
for a = 1:length(xss)-2;
    rot(a) = abs(angle(a)-angle(a+1));
    dist(a) = sqrt((xss(a)-xss(a+2)).^2 + (yss(a)-yss(a+2)).^2);
end
%     dist = dist(200:end-199);
rot(rot > 180) = rot(rot > 180)-180;
rot = 360-rot; %want rotation to be small so fixation values are all small
%     rot = rot(200:end-199);

points = [dist' vel' accel' rot'];
for ii = 1:size(points,2) %normalizes points to [0 1] by parameter
%     thresh = nanmean(points(:,ii))+nanstd(points(:,ii));%move outliers
%     points((points(:,ii) > thresh),ii) = thresh;
    points(:,ii) = points(:,ii)-min(points(:,ii));
    points(:,ii) = points(:,ii)/max(points(:,ii));
end

[r,c]=find(isnan(points));
points_nonan = points(setxor(1:size(points,1),unique(r)),:);
    
cd('R:\Mike\MATLAB\ClusterFix')

% %---Global Clustering---%
% sil = zeros(1,5); %determines the number of clusters by comparing the ratio
% %of intercluster and intracluster distances, faster mod of silhouette
% for numclusts = 2:5
%     fprintf('Now analyzing %g clusters\n',numclusts)
%     T = kmeans(points_nonan(1:10:end,2:4),numclusts,'replicate',5);
%     [silh] = InterVSIntraDist(points_nonan(1:10:end,2:4),T);
%     sil(numclusts) = mean(silh);
% end
% sil(sil > 0.9*max(sil)) = 1;
% numclusters = find(sil == max(sil));
% T = kmeans(points_nonan,numclusters(end),'replicate',10);

% just use 5 clusters
T = kmeans(points_nonan,5,'replicate',10);
Torig = T;
    
meanvalues = zeros(max(T),size(points_nonan,2));
stdvalues = zeros(max(T),size(points_nonan,2));
for TIND = 1:max(T)
    tc = find(T == TIND);
    meanvalues(TIND,:) = mean(points_nonan(tc,:));
    stdvalues(TIND,:) = std(points_nonan(tc,:));
end

% determines fixation clusters by overlapping distributions in velocity
% and acceleration state space, here assumes gaussian distributions
[~, fixationcluster] = min(sum(meanvalues(:,2:3),2));
T(T == fixationcluster) = 100;
fixationcluster2 = find(meanvalues(:,2) < meanvalues(fixationcluster,2)...
    +3*stdvalues(fixationcluster,2));
fixationcluster2(fixationcluster2 == fixationcluster)= [];
for iii = 1:length(fixationcluster2)
    T(T == fixationcluster2(iii)) = 100;
end
T(T ~= 100) = 2;
T(T == 100) = 1;

figure;hold on
plot(points_nonan(:,2))
scatter(1:size(points_nonan,1),T/2,'or')
% xlim([157565 162772])
    
fixationindexes =  find(T == 1)';
[fixationtimes] = BehavioralIndex(fixationindexes);
fixationtimes(:,(diff(fixationtimes,1) < 25)) = []; %25 ms duration threshold

figure
ax(1)=subplot(2,1,1);
hold on
plot(points_nonan(1:10000,2))
for k=1:find(fixationtimes(2,:)<=10000,1,'last')
    line([fixationtimes(1,k) fixationtimes(1,k)],ylim,'Color','g')
    line([fixationtimes(2,k) fixationtimes(2,k)],ylim,'Color','r')
end
ax(2)=subplot(2,1,2);
h=[xss(1:20000);yss(1:20000)];
h=h(:,~isnan(h(1,:)));
plot(1:10000,h(:,1:10000))
linkaxes(ax,'x')
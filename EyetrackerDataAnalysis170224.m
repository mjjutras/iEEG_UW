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



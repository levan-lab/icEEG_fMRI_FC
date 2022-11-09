% computes and plots static BOLD FC
clear all
close all

subject = 'sub-049';
run = '02';
prepro_dir = ['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject];
ts_dir = [prepro_dir,'/rs_fc/ts/24hmp+2p_bpbutter/'];
load([prepro_dir,'/conn_analysis/',subject,'_electrodes_rois_eegfmri_select.mat'])

numRois = height(rois)

cd(ts_dir)

file = {};
roi_pairs = {};
regions = {};
hemi = {};
regions_hemi = {};
ts = [];

for i = 1:numRois
    file = [rois.rois{i},'.txt']
    roi_pairs{i,1} = rois.rois(i);
    regions = cellstr(rois.regions);
    hemi = cellstr(rois.hemi);
    regions_hemi{i,1} = [regions{i,1},'_',hemi{i,1}];
%    rois_name(i,1) = rois(i,1(1:end-4);
    ts(:,i) = dlmread(file);

end

contacts_corr = corr(ts);

% Colormap 1 (hot)
T = [1 0 0   % Red   
    1 1 0]; % Yellow                 
x = [0 1];
x = x(end:-1:1);
c1 = interp1(x,T,linspace(0,1,64));

% Colormap 2 (cold)
T = [0 1 1             % Turqoise
     0 0.0745 0.6078]; % Dark blue         
x = [0 1];
x = x(end:-1:1);
c2 = interp1(x,T,linspace(0,1,64));        

% Combine and add white
c = [c2;c1];
c(64:65,:) = 1;


x1= 50;% position (in pixels) of the lower left corner on your screen in x direction
y1= 50; % position (in pixels) of the lower left corner on your screen in y direction
xsize= 1700; % size of your figure in x-direction
ysize= 700; % size of your figure in y-direction

figure('Position',[x1,y1,xsize,ysize])

subplot(1,2,1)
% Plot
imagesc(contacts_corr,[-1 1]);
colorbar;
colormap(c)


axis = regions_hemi;

set(gca, 'XTick', [1:1:numRois], 'XTickLabel', axis, 'YTick', [1:1:numRois], 'YTickLabel', axis)

xtickangle(90)

title(['Correlations ', subject,'; run ', run,' 24HMP+2P bpf'])

subplot(1,2,2)
% Plot
imagesc(contacts_corr,[-1 1]);
colorbar;
colormap(c)


axis = roi_pairs;

set(gca, 'XTick', [1:1:numRois], 'XTickLabel', axis, 'YTick', [1:1:numRois], 'YTickLabel', axis)

xtickangle(90)

title(['Correlations ', subject,'; run ', run,' 24HMP+2P bpf'])

%%%%

B = ones(length(contacts_corr),length(contacts_corr));
tri= triu(contacts_corr);
bold_fc = tri(logical(triu(B,1)));

%%%%%
cd ([prepro_dir,'/conn_analysis/'])

save(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject,'/conn_analysis/run',run,'_24hmp2p_butter_rois_fc_select.mat'],'ts','contacts_corr','roi_pairs','regions_hemi','bold_fc')


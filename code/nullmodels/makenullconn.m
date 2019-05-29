cd(homedir); addpath(genpath(pwd));
load(fullfile(savedir,'W.mat'));

Wnull = randmio_dir(W,10);

save(fullfile(savedir,'Wnull.mat'),'Wnull');

coor = csvread('Data83018/coor.csv',1,1);
D = squareform(pdist(coor,'Euclidean'));    % R first then L
save(fullfile(savedir,'D.mat'),'D')

%% degree sequence, weight-length, weight, length preserving null
nbins = 5; nrewire=1e5; 
%[~,Wrw] = fcn_preserve_degseq_lengthdist(W,D,nbins,nrewire);

%% degree sequence preserving rewiring
Wrw = dir_generate_srand(W,1e6);

save(fullfile(savedir,'Wrewire.mat'),'Wrw');

%check degree is same
%plot(sum(W,2),sum(Wrw,2),'.')

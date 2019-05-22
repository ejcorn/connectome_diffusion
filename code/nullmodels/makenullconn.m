clear all; close all; clc
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/BCT/'));
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/brainmapping2/Colormaps'));
addpath(genpath('~/Dropbox/Neurodegeneration/MouseDiffusion/'));
cd ~/Dropbox/Neurodegeneration/MouseDiffusion/
load([homedir,'processed/W.mat']);

Wnull = randmio_dir(W,10);

save([savedir,'Wnull.mat'],'Wnull');

coor = csvread('Data83018/coor.csv',1,1);
D = squareform(pdist(coor,'Euclidean'));    % R first then L
save([savedir,'D.mat'],'D')

%% degree sequence, weight-length, weight, length preserving null
nbins = 5; nrewire=1e5; 
%[~,Wrw] = fcn_preserve_degseq_lengthdist(W,D,nbins,nrewire);

%% degree sequence preserving rewiring
Wrw = dir_generate_srand(W,1e6);

save([savedir,'Wrewire.mat'],'Wrw');

plot(sum(W,2),sum(Wrw,2),'.')

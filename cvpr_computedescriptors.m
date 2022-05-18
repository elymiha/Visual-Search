%% EEE3032 - Computer Vision and Pattern Recognition (ee3.cvpr)
%%
%% cvpr_computedescriptors.m
%%
%% This code will iterate through every image in the MSRCv2 dataset
%% and call a function 'extractDescriptor' to extract a descriptor from the
%% image. 
%%
%% (c) John Collomosse 2010  (J.Collomosse@surrey.ac.uk)
%% Centre for Vision Speech and Signal Processing (CVSSP)
%% University of Surrey, United Kingdom
%% (c) Elena Mihalas 2016 (em00561@surrey.ac.uk)

close all;
clear all;

%% Folder that holds the collection of images
DATASET_FOLDER = 'MSRC_ObjCategImageDatabase_v2';

%% Create a folder to hold the results...
OUT_FOLDER = 'C:/Users/Elena/Desktop/Università/Computer Vision and Pattern Recognition/Lab/Matlab codes/Assignment/descriptors';
%% and within that folder, create another folder to hold these descriptors
%% the idea is all your descriptors are in individual folders - within
%% the folder specified as 'OUT_FOLDER'.

OUT_SUBFOLDER='globalRGBhisto';
F_all=[];
allfiles=dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));

for filenum=1:length(allfiles)
    fname=allfiles(filenum).name;
    % tic;
    imgfname_full=([DATASET_FOLDER,'/Images/',fname]);
    img = imread(imgfname_full);
    fout=[OUT_FOLDER,'/',OUT_SUBFOLDER,'/',fname(1:end-4),'.mat']; %replace .bmp with .mat
    [F]=extractDescriptor(img);
    save(fout,'F');
    % toc
    F_all=[F_all; F];
   
end

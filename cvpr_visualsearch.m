%% EEE3032 - Computer Vision and Pattern Recognition (ee3.cvpr)
%%
%% cvpr_visualsearch.m
%%
%% This code will load in all descriptors pre-computed (by the
%% function cvpr_computedescriptors) from the images in the MSRCv2 dataset.
%%
%% It will pick a particular descriptor and compare all other descriptors to
%% it - by calling cvpr_compare.  In doing so it will rank the images by
%% similarity to the selected descriptor.  Note that by default the
%% function cvpr_compare returns the Euclidean distance, but it can be 
%% programmed to return some other distance metric such as L1 norm and
%% mahalanobis distance.
%%
%% (c) John Collomosse 2010  (J.Collomosse@surrey.ac.uk)
%% Centre for Vision Speech and Signal Processing (CVSSP)
%% University of Surrey, United Kingdom
%% (c) Elena Mihalas 2016 (em00561@surrey.ac.uk)

close all;
clear all;

%% Folder that holds the collection of images
DATASET_FOLDER = 'MSRC_ObjCategImageDatabase_v2';

%% Folder that holds the results...
DESCRIPTOR_FOLDER = 'descriptors';
%% and within that folder, another folder to hold the descriptors
%% we are interested in working with
DESCRIPTOR_SUBFOLDER='globalRGBhisto';


%% 1) Load all the descriptors into "ALLFEAT"
%% each row of ALLFEAT is a descriptor (is an image)

ALLFEAT  = [];
ALLFILES = cell(1,0);
ctr      = 1;
allfiles = dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));

for filenum = 1:length(allfiles)
    fname         = allfiles(filenum).name;
    imgfname_full =([DATASET_FOLDER,'/Images/',fname]);
    img           = double(imread(imgfname_full))./255;
    thesefeat     = [];
    featfile      = [DESCRIPTOR_FOLDER,'/',DESCRIPTOR_SUBFOLDER,'/',fname(1:end-4),'.mat']; %replace .bmp with .mat
    load(featfile,'F');      % Loads the color histograms
    ALLFILES{ctr} = imgfname_full;
    ALLFEAT       = [ALLFEAT ; F];   % Puts in a column all the histograms
    ctr           = ctr+1;
end




%% Choose between a _single_ query search or a _complete_ search with respect to all possible queries 

NIMG        = size(ALLFEAT,1);  % number of images in collection
top         = 15;               % Top n results performance
dim         = 15;               % Dimension of the descriptor matrix
type_search = 'complete';


switch (type_search) 
    
%% 2) Compute the distance of image to the query + Precision Recall curve    
    case 'single'               % A single_query search
     
    % Pick a particular image from 591 images to be the query   
    queryimg   = 521 ;
    dst        = [];
    % Apply the Principal Component Analysis
    e          = Eigen_Build(ALLFEAT');
    
    % Project the image descriptor into a lower dimentional space
     ALLFEATPCA = descriptor_projection( ALLFEAT', e, dim ); % I obtain a [3x591] matrix
     ALLFEATPCA = ALLFEATPCA';
  
    % Compute the distance of the query image from each image
    for i=1:NIMG
        candidate = ALLFEATPCA(i,:);
        query     = ALLFEATPCA(queryimg,:);
        thedst    = cvpr_compare(query,candidate, e.val, dim); % Query and candidate are two histograms
        dst       = [dst ; [thedst i]];
    end
    
    % Compute the precision-recall curve of single query images upon "n"
    % top closest distances
    dst      = sortrows(dst,1);                                % sort the results in ascending order for the first column
    [R,P,AP] = Precision_recall(queryimg, ALLFILES, dst, top, type_search);
    
    namer = ['Recall ' 'top' num2str(top)];
    namep = ['Precision ' 'top' num2str(top)];
    
    % Visualise the results
    f1 = figure;
    plot(R,P, 'r')
    hold on
    plot(R,P,'r*')
    xlabel(namer)
    ylabel(namep)
    title('Precision-Recall Curve')
    legend(['AP=' num2str(AP)])
    
    % Save the result
    name1 = ['HARRIS_PR_' num2str(queryimg) '_top' num2str(top) '_city_block' '.png'];
    saveas(f1, name1)
    
    % Show the first 10 closest images to the query image
    SHOW       = 10;             
    dst        = dst(1:SHOW,:);
    outdisplay = [];
    
    for i=1:size(dst,1)
        img        = double(imread(ALLFILES{dst(i,2)}))./255;
        img        = img(1:2:end,1:2:end,:); % make image a quarter size
        img        = img(1:95,:,:); % crop image to uniform size vertically 
        outdisplay = [outdisplay img];
    end
    
    f2 = figure;
    imshow(outdisplay);
    axis off;
    
    % Save the picture
    name2 = ['HARRIS_img_' num2str(queryimg) '_city_block' '.png'];
    saveas(f2, name2)
    
    
%% 3) Compute the distance of image to the query + MAP + Confusion Matrix
    case 'complete'            % a complete search with respect to all possible queries
        
        AP_all    = [];
        Precision = [];
        Recall    = [];
        
        % Number of Categories (20)

s_obj = {'meadow', 'tree', 'house', 'airplane', 'cow', 'face', 'car', ...
         'bicycle', 'sheep', 'flower', 'signboard', 'bird', 'books'...
         'chair', 'cat', 'dog', 'street', 'water', 'people', 'sea'};  
     
       
        n_objects   = size(s_obj,2);
        mat_conf  = zeros(n_objects, n_objects);
        
        for j=1:NIMG
            queryimg   = j;
            dst        = [];
            % Apply the Principal Component Analysis
            e          = Eigen_Build(ALLFEAT');
            
            % Project the image descriptor into a lower dimentional space
            % in order to obtain a [dim x 591] matrix
            ALLFEATPCA = descriptor_projection( ALLFEAT', e, dim ); 
            ALLFEATPCA = ALLFEATPCA';
           
            % Compute the distance of the query image from each image
            for i=1:NIMG
                candidate = ALLFEATPCA(i,:);
                query     = ALLFEATPCA(queryimg,:);
                thedst    = cvpr_compare(query,candidate, e.val, dim);   
                dst       = [dst ; [thedst i]];
            end
            
            dst        = sortrows(dst,1);   % sort the results in ascending 
                                            % order for the first column
            
            % Compute the Average Precision for each query image
            [R,P,AP]   = Precision_recall(queryimg, ALLFILES, dst, top, type_search);
            AP_all     = [AP_all, AP];
            Precision  = [Precision, mean(P)];
            Recall     = [Recall, mean(R)];
            
            % Compute the Confusion Matrix
            mat_conf = Confusion_matrix(queryimg, ALLFILES, mat_conf, dst, top);
            
        end
        
        % Plot the Confusion Matrix
        norm_conf = sum(mat_conf,2);          % Sum the values in each row
        mat_norm  = repmat(norm_conf,1,20);   
        mat_conf  = mat_conf./mat_norm;       % Normalisation step - divide 
                                              % each element in the row by 
                                              % the sum across that row
        
        MAP_mat_conf = sum(diag(mat_conf))/20 % Mean Average Precision of the 15 top distances
    
        f3 = figure;
        imagesc(mat_conf)     
        %colormap(flipud(gray));
    
        title(['Confusion Matrix for the top ' num2str(top) ' distance values'])
        axis image
        set(gca, 'XTick', [1,5,8,12,16,20]);
        set(gca, 'XTickLabel', {'meadow', 'cow', 'bicycle', 'bird', 'dog', 'sea'});
        set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]);
        set(gca, 'YTickLabel', s_obj);
        set(gca, 'FontSize', 8.0);
        xlabel('Intended Label')
        ylabel('Ground Truth Label')
        colorbar    
        
        % Save the figure
        saveas(f3, 'GCH_Confusion Matrix_DIM15_mah_q7.png')
        
        % Compute the Mean Average Precision upon all the images in the
        % dataset
        MAP       = mean(AP_all);
        Recall    = sort(Recall, 2, 'ascend');
        Precision = sort(Precision, 2, 'descend');
        
        % Make the values inside Recall and Precision vary from 0 to 1
        Min_R     = repmat(min(Recall), 1, length(Recall));
        Min_P     = repmat(min(Precision), 1, length(Precision));
        Recall    = (Recall-Min_R)./(max(Recall)-min(Recall));
        Precision = (Precision-Min_P)./(max(Precision)-min(Precision));
        
        % Visualise the results
        f4 = figure;
        plot(Recall,Precision, 'r')
        xlabel('Average Recall')
        ylabel('Average Precision')
        title('Precision-Recall Curve for the Harris Corner Descriptor')
        legend(['MAP=' num2str(MAP)])
       
        % Save the result
        saveas(f4, 'GCH_MAP_DIM15_mah_q7.png')

        
end
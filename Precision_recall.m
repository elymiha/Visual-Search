function [R, P, AP] = Precision_recall(query, ALLFILES, dst, top, type_search )

% This function computes the precision and recall rate of individual query
% image. Moreover it evaluates the average precision.

N           = length(ALLFILES);     % Number of images in the dataset
GroundTruth = zeros(N,1);           % Size of the ground truth [591x1]
P           = [];
R           = [];

q = ALLFILES{query};    % Select the query image from the dataset
q = q(140:(end-6));     % Extract the number of the image

for n=1:N
    each_img = ALLFILES{n};
    each_img = each_img(140:(end-6));
    
    if  (each_img(1:2) == q(1:2))    % The first two characters identify the category of the image. 
        GroundTruth(n,1) = 1;        % Create the mask
    end
    
end

relevant  = nnz(GroundTruth);        % Equal with number of elements in the category to which the query image belongs
sorted_GT = GroundTruth(dst(:,2));   % Sorts the elements in Ground Truth with respect to the 2nd column of the "dst" matrix
count     = 0;

% Differentiate between complete search and single search 
% in order to specify the top "n" elements for the single case

switch(type_search)
    case 'single'
        for j=2:34                   % The double loop makes count=594
            for i=1:j
                count = count+1;
                
                if count>top                % The count should not exceed the number of files in the dataset
                    break                   % It is possible to define count as top 15 evaluation
                end
                
                returned_rel = nnz(sorted_GT(1:count,1));  % Number of relevant returned elements within progressively increasing
                P            = [P; returned_rel/count];    % Precision computation
                R            = [R; returned_rel/relevant]; % Recall computation
                
            end
        end
        
        
        AP = sum((P.* sorted_GT(1:top,1))/relevant);                % Average Precision computation
        
    case 'complete'
        for j=2:34                          % The double loop makes count=594
            for i=1:j
                count = count+1;
                
                if count>591                % The count should not exceed the number of files in the dataset
                    break                  
                end
                
                returned_rel = nnz(sorted_GT(1:count,1));  % Number of relevant returned elements within progressively increasing
                P            = [P; returned_rel/count];    % Precision computation
                R            = [R; returned_rel/relevant]; % Recall computation
                
            end
        end
        
        
        AP = sum((P.* sorted_GT)/relevant);                % Average Precision computation
        
end
end


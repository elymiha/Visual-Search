function [ mat_conf ] = Confusion_matrix( queryimg, ALLFILES,mat_conf, dst,top )

% This function computes the confusion matrix for the top 15 results

 % Number corresponding to each category 
classes = {'1_', '2_', '3_', '4_', '5_', '6_', '7_', '8_', '9_', '10',...
           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'};
       
n_cat =  size(classes,2);            
query = ALLFILES{queryimg};  % Select the query image from the dataset
query = query(140:(end-6));  % Extract the number of the image
q     = ismember(classes, query(1:2));

if any(q)
    main_category = find(q); % Discover to which category the query image belongs
end


for j=1:top
    each_img = ALLFILES{dst(j,2)};                % Find the category of the first 15 top closest images
    each_img = each_img(140:(end-6));
    member   = ismember(classes, each_img(1:2));  % The first two characters identify the category of the image.
    
    if any(member)
        category = find(member);
        mat_conf(main_category,category)= mat_conf(main_category, category)+1;
    end
end

end


function dst=cvpr_compare(F1, F2, variance, dim)

% This function compare F1 to F2 by computing the distance between the two descriptors
% To switch between the computation of different distances just assign one of
% the options (euclidean, city_block, mahalanobis) to distance

distance = 'mahalanobis';
 
switch ( distance )
    case 'euclidean'
        
    dst = sqrt(sum((F1-F2).^2));   % d= sqrt(f_1.^2+ f_2.^2 )
    
    case 'city_block'
        
    dst = sum(abs(F1-F2));         % d= f1 + f2
    
    case 'mahalanobis'
        
    % Note that when you select this distance metric, you should uncomment
    % the ALLFEATPCA piece of code in "cvpr_visualsearch" and you should
    % specify query and candidate as a search through ALLFEATPCA and NOT
    % ALLFEAT.
    
   
    variance(variance==0)=1;       % Check for eigenvalues of 0 and equal them to one
    
    variance = variance(1:dim);
    dst   = sqrt(sum((F1-F2).*(F1-F2)./variance'));
    
    
    
end
        
  
        
    


return;

function [F] = extractDescriptor(img)

% Returns a row [...] representing an image descriptor computed from image 'img'

% To switch between the computation of different descriptors just assign one of
% the options (global_colour, spacial_colour, texture, colour and texture, harris) to descriptor

descriptor = 'global_colour' ;

switch(descriptor)
    
    case 'global_colour'
        Q    = 7;
        qimg = double(img)./256;    % Normalised RGB image in range 0 to (Q-1)
        qimg = floor(qimg.*Q);      % Drop the decimal point
        
        % Pixel's RGB value summarised by a single integer value in range [0:Q^3-1]
        
        bin  = qimg(:,:,1)*Q^2 + qimg(:,:,2)*Q + qimg(:,:,3); % 2D image
        
        % Build a frequency histogram from this values
        
        vals = reshape(bin,1,size(bin,1)*size(bin,2)); % obtain a long vector
        
        H    = hist(vals, Q^3);
        
        F    = H./sum(H);                              % Normalisation of H
 
%_____________________________________________________________________
    case 'spacial_colour'
        
        gimg               = double(img)./255;
        [rows, columns,c ] = size(gimg);
        a                  = floor(rows/6);       % Choose 6x6 macroblocks 
        b                  = floor(columns/6);
        F                  = [];
        
        d2 = 0;
        for j=1:6       % Get along the columns
            
            d1=0;
            for i=1:6   % Get along the rows
                
            szw   = [1:a]+d1*a;   % d1 changes inside the i loop
            szl   = [1:b]+d2*b;   % d2 changes inside the j loop
            image = gimg(szw,szl,:);
            
            R = mean(mean(image(:,:,1)));
            G = mean(mean(image(:,:,2)));
            B = mean(mean(image(:,:,3)));
            F = [F, [R,G,B]];
            
            d1 = d1+1;  
            end
            
            d2 = d2+1;
        end
 %____________________________________________________________________       
    case 'texture'
        
        gimg    = double(img)./255;
        greyimg = gimg(:,:,1)*0.30 + gimg(:,:,2)*0.59 + gimg(:,:,3)*0.11; 
        
        [rows,columns] = size(greyimg);
        a              = floor(rows/6);
        b              = floor(columns/6);
        F              = [];
        d2             = 0;
        
        % This double loop iterates vertically down along the columns
        for j=1:6       % Get along the columns
            
            d1 = 0;
            for i=1:6   % Get along the rows
                
            szw   = [1:a]+d1*a;   % d1 changes inside the i loop
            szl   = [1:b]+d2*b;   % d2 changes inside the j loop
            image = greyimg(szw,szl);
            
            sobel = [1 0 -1; 2 0 -2; 1 0 -1];
            
            % compute horizontal and vertical gradients
            Ix = conv2(image, sobel, 'same');
            Iy = conv2(image, sobel', 'same');
            
            theta       = atan2(Iy,Ix);          % Gives numbers in the range [-pi:pi]
            mag         = sqrt(Ix.^2 + Iy.^2);
            thresholded = double((mag > 0.25));  % I ignore all the data that have a magnitude less than 0.25
            
            relevant_theta = theta.*thresholded; 
            vals           = nonzeros(relevant_theta);
          
            H            = hist(vals, 8);        % bins the elements of vals into 8 equally spaced containers
            Normalised_H = H./sum(H);
            F            = [F, Normalised_H];
          
            d1 = d1+1;
            end
            
            d2 = d2+1;
        end
        
        % Good practice to check the presence of eventual NaN values
        if (nnz(isnan(F))>0)
            idx   = find(isnan(F));
            F(idx)= 0;
        end
        
%_______________________________________________________________________        
    case 'colour_and_texture'
        gimg    = double(img)./255;
        greyimg = gimg(:,:,1)*0.30 + gimg(:,:,2)*0.59 + gimg(:,:,3)*0.11; 
        
        [rows, columns, c ] = size(gimg);
        a                   = floor(rows/6);
        b                   = floor(columns/6);
        F                   = [];
        
        d2 = 0;
        for j=1:6       % Get along the columns
            
            d1 = 0;
            for i=1:6   % Get along the rows
                
                szw   = [1:a]+d1*a;   % d1 changes inside the i loop
                szl   = [1:b]+d2*b;   % d2 changes inside the j loop
                image = gimg(szw,szl,:);
                
                R = mean(mean(image(:,:,1)));
                G = mean(mean(image(:,:,2)));
                B = mean(mean(image(:,:,3)));
                
                image1 = greyimg(szw,szl);
                sobel  = [1 0 -1; 2 0 -2; 1 0 -1];
                
                % compute horizontal and vertical gradients
                Ix = conv2(image1, sobel, 'same');
                Iy = conv2(image1, sobel', 'same');
                
                theta       = atan2(Iy,Ix);
                mag         = sqrt(Ix.^2 + Iy.^2);
                thresholded = double((mag >0.25));   % I ignore all the data that have a magnitude less than 0.25
                
                relevant_theta = theta.*thresholded; 
                vals           = nonzeros(relevant_theta);
                
                H            = hist(vals, 8);        % bins the elements of vals into 8 equally spaced containers
                Normalised_H = H./sum(H);
                F            = [F, Normalised_H,[R,G,B]];
                
                d1 = d1+1;
            end
            
            d2 = d2+1;
        end
        
        % Good practice to check the presence of eventual NaN values
        if (nnz(isnan(F))>0)
            idx   = find(isnan(F));
            F(idx)= 0;
        end
        
%________________________________________________________________________              
    case 'harris'
        
        gimg    = double(img)./255;
        
        % Run the Harris corner detector
        thresh=200; % top 1000 corners
        F = harris(gimg, thresh);
        
       
end

return;
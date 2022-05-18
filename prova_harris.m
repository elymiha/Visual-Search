
close all
clear all

img=imread('/MSRC_ObjCategImageDatabase_v2/Images/13_12_s.bmp');
gimg= double(img)./255;
greyimg = gimg(:,:,1)*0.30 + gimg(:,:,2)*0.59 + gimg(:,:,3)*0.11; 

% Run the Harris corner detector
ncorners=200; 
width = 3; 
sigma = 1;
subpixel = 0;
 
mask = [-1 0 1; -1 0 1; -1 0 1] / 3;

% compute horizontal and vertical gradients
Ix = conv2(greyimg, mask, 'valid');
Iy = conv2(greyimg, mask', 'valid');

% These elements encode local intensity variation
Ixy = Ix .* Iy;
Ix = Ix.^2;
Iy = Iy.^2;

% smooth them
gmask = torr_gauss_mask(width, sigma);

Ix = conv2(Ix, gmask, 'valid');
Iy = conv2(Iy, gmask, 'valid');
Ixy = conv2(Ixy, gmask, 'valid');

% Compute corner strength (cornerness)
c = (Ix + Iy) - 0.04 * (Ix .* Iy - Ixy.^2);
[r,col] = size(c);

% Select 3x3 windows inside the image and compute the max 
c_pad = padarray(c, [1,1]);
[rows,columns] = size(c_pad);
nbr = floor(rows/3);
nbc = floor(columns/3);
w_r = floor(rows/nbr);
w_c = floor(columns/nbc);

d2      = 0;
mat_max = zeros(r, col);

        % This double loop iterates vertically down along the columns
        for j=1:r       % Get along the columns
            
            d1 = 0;
            for i=1:col   % Get along the rows
                
                szw   = [1:w_r]+d1;   % d1 changes inside the i loop
                szl   = [1:w_c]+d2;   % d2 changes inside the j loop
                if (szw(1,end)>columns || szl(1,end)>rows)  
                    break
                end
                
                image = c_pad(szl,szw);
                mat_max(j,i) = max(max(image));
                
                d1 = d1+1;
                
            end
            d2 = d2+1;
        end
        
% if pixel equals max, it is a local max, find index
ci3 = find(c == mat_max); % column of max indexes
cs3 = c(ci3);             % column of max values

[cs2,ci2] = sort(cs3); %ci2 2 is an index into ci3 which is an index into c


% put strongest ncorners in a list cs together with indices ci

l = length(cs2);
cs = cs2(l-ncorners+1:l);
ci2s = ci2(l-ncorners+1:l);

ci = ci3(ci2s);

corn_thresh = cs(1);

% row and column of each corner
[nrows, ncols] = size(c);


c_row = rem(ci,nrows) +(width+1);
c_col = ( ci - c_row  +width+1)/nrows + 1 +width+1;

% to convert to x,y we need to convert from rows to y
c_coord = [c_col c_row];

if subpixel

    %   Extension to subpixel accuracy - John Collomosse 2002 jpc@cs.bath.ac.uk
    
    %   for each 'good corner' we must refine to subpixel accuracy, we fit
    %   a quadratic surface and find zero grad of that surface
           
    % 1) Fit the quadratic surface to a 3x3 neighbourhood centered around the 'corner' pixel
    %    (least squares fit)
    for lp=1:size(c_coord,1)
        x=c_coord(lp,1)-(width+1);    % convert into index into c (which is smaller than image)
        y=c_coord(lp,2)-(width+1); 
        
        if (x==0 || x==1)
            x=2;
        end
        
        if (y==0 || y==1)
            y=2;
        end
        
        A=[x^2      y^2       x*y           x     y     1   -c(y,x); ...
          (x-1)^2   (y-1)^2   (x-1)*(y-1)   x-1   y-1   1   -c(y-1,x-1); ...
          (x+1)^2   (y-1)^2   (x+1)*(y-1)   x+1   y-1   1   -c(y-1,x+1); ...
          (x-1)^2   (y+1)^2   (x-1)*(y+1)   x-1   y+1   1   -c(y+1,x-1); ...
          (x+1)^2   (y+1)^2   (x+1)*(y+1)   x+1   y+1   1   -c(y+1,x+1); ...
          x^2       (y+1)^2   x*(y+1)       x     y+1   1   -c(y+1,x); ...
          x^2       (y-1)^2   x*(y-1)       x     y-1   1   -c(y-1,x); ...
          (x-1)^2   y^2       (x-1)*y       (x-1) y   1     -c(y,x-1); ...
          (x+1)^2   y^2       (x+1)*y       (x+1) y   1     -c(y,x+1)];
        [u s v]=svd(A);
        res=v(:,end);
        res=res./res(end,1);

        % biquad. surface fitted, here are the params
        % I(x,y)=ax^2+by^2+cxy+dx+ey+f
        a=res(1);
        b=res(2);
        c2=res(3); % note c2, we used 'c' for something else
        d=res(4);
        e=res(5);
        f=res(6);
        
        coord=inv([2*a c2 ; c2 2*b])*[-d ; -e];        % solve to find the turning pt of surface
        
        % in some, rare, cases the max/min will be outside of the 3x3 block we fitted the
        % surface to. In this case, revert to single pixel accuracy.. (unclear to me why
        % this happens).
        if abs(x-coord(1))>1
            coord(1)=x;
        end
        if abs(y-coord(2))>1
            coord(2)=y;
        end
    
        coord=coord+width+1;  % corners converted back to image space
        c_coord(lp,1:2)=coord';
         
    end
end

% Extract the values from a 3x3 window fitted on each corner
% The matrices X and Y are used to select the values local to a corner without
% modifing them
X = [0 0 0; 0 1 0; 0 0 0];   
Y = [ 0 0 0; 0 1 0; 0 0 0];
F =[];                      % It will hold all mean R,G,B values local to the corners

for i=1:size(c_coord,1)     % I consider all top "n" corners
    RGB = [0,0,0];
for j=1:3
    c_row       = c_coord(i,2);   % Select the row
    c_col       = c_coord(i,1);   % Select the column
    c_col_left  = c_col-1;        % Select the adjacent left column
    c_col_right = c_col+1;        % Select the adjacent right column
    
    % Each of the Z(i) extract a column of three values local to the corner
    Z1 = gimg(c_row-1:1:c_row+1,c_col_left-1:1:c_col_left+1,j)*X*Y;
    Z2 = gimg(c_row-1:1:c_row+1,c_col-1:1:c_col+1,j)*X*Y;
    Z3 = gimg(c_row-1:1:c_row+1,c_col_right-1:1:c_col_right+1,j)*X*Y;
    % Concatenate the three [3x1] columns in order to obtain the 3x3 window
    % with image colour values local to the corner
    Z  = [Z1(:,2),Z2(:,2),Z3(:,2)];
    
    % Take the mean for each of the three colours
    RGB(j) = mean(mean(Z(:,:)));  
end
    F = [F, RGB];  
    
end

% Plot positions of the Harris corners
figure;
imshow(gimg);
hold on;
plot(c_coord(:,1),c_coord(:,2),'rx');



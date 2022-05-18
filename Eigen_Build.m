%% EEE3032 - Computer Vision and Pattern Recognition (ee3.cvpr)
%%
%% Eigen_Build.m
%% Builds an eigen model from n observation of d-dimensional data  
%%
%% Usage:  E = Eigen_Build(observations)
%%
%% IN:  observations   - d row, n column matrix of observations
%%
%% OUT: E              - Eigenmodel structure
%%                       E.org - mean
%%                       E.vct - matrix of eigenvectors one per column 
%%                       E.val - column of eigenvalues(matching E.vct cols)
%%                       E.N   - number of observations used i.e. = n
%%
%% (c) John Collomosse 2010  (J.Collomosse@surrey.ac.uk)
%% Centre for Vision Speech and Signal Processing (CVSSP)
%% University of Surrey, United Kingdom

function E=Eigen_Build(obs)

    E.N  =size(obs,2);  % takes the number of columns or n -> 591
    E.D  =size(obs,1);  % takes the number of dimensions -> histogram length
    E.org=mean(obs')';  % takes the mean along each dimension
    
    % Computation of the std deviation
    obs_translated=obs-repmat(E.org,1,E.N); % Repmat creates a [d x n] matrix
                                            % containing the mean values
                                            % along each direction, then I
                                            % substract from each original
                                            % pixel the corresponding mean
    
    % Covariance matrix
    C=(1/E.N) * (obs_translated * obs_translated'); 
    
    % Eigenvalue Decomposition
    [U V]=eig(C);   % U= matrix of eigenvectors
                    % V= matrix of eigenvalues
    
    % sort eigenvectors and eigenvalues by eigenvalue size (desc)
    linV=V*ones(size(V,2),1);  % Obtain a [n x 1] matrix -> sum the values along the row (->)
    S=[linV U'];               % Obtain a [n x d+1] matrix
    S=flipud(sortrows(S,1));   % SORTROWS(X,COL) sorts the matrix based on 
                               % the columns specified in the vector COL.
                               % FLIPUD(X) returns X with columns preserved
                               % and rows flipped in the up/down direction
    U=S(:,2:end)';
    V=S(:,1);
    
    E.vct=U;
    E.val=V;
    
function [ ALLFEATPCA] = descriptor_projection( descriptors, e, dim)
% This function plots the descriptor into a lower dimensional space

% Estimate how much of the total energy is present in the first [1:dim]
% eigenvalues
totalenergy    = sum(abs(e.val));
current_energy = sum(abs(e.val(1:dim)))/totalenergy;
current_energy = round(current_energy*100);

%fprintf('To %d dimensions corresponds circa %d percent of the total energy\n',dim, current_energy);

e.val     = e.val(1:dim);
e.vct     = e.vct(:,1:dim);

descriptors_translated = descriptors-repmat(e.org,1,size(descriptors,2));
ALLFEATPCA             = e.vct'*descriptors_translated;


% If dim=3 it is possible to visualise the descriptors' distribution 

% figure;
% plot3(ALLFEATPCA(1,:), ALLFEATPCA(2,:), ALLFEATPCA(3,:), 'b*');
% xlabel('1st eigenvector');
% ylabel('2nd eigenvector');
% zlabel('3rd eigenvector');

end


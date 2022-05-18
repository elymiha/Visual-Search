# Visual-Search
A Matlab visual search program, which evaluates the implementation of different approaches for images similarity estimation and different distance measures

The visual search matlab code is composed of two main scripts, namely the cvpr_computedescriptors 
and cvpr_visualsearch.

1) cvpr_computedescriptors -> it iterates through each image in the MSRC_... folder and creates
an image descriptor by calling the function "extractDescriptor.m". The latter function generates the
specified descriptor by computing an histogram of a specific content which is present inside each image.
In fact it allows to generate 5 different histograms:
- global  colour histogram
- spacial colour distribution histogram
- spacial edge orientation histogram
- spacial colour and edge orientation histogram
- harris  corner descriptor -> "harris" function

2) cvpr_visualsearch -> it loads all the descriptors computed by "cvpr_computedescriptors" into ALLFEAT matrix
Then the script is divided into 2 parts. The first part implements the visual search for a particular query image.
Thus it computes, through the function "cvpr_compare"  the distance between the query descriptor and the 
descriptor from each image. Note that it is possible to select between three different distances 
(city block, euclidean, mahalanobis ).Furthermore it evaluates how discriminative and robust is the selected
descriptor by applying the "Precision_recall.m" function.
The second part, instead, implements a complete search by considering for once each image in the dataset a query image. 
It also computes the distance between 2 particular descriptors at a time and it evaluates the similarity by computing
the confusion matrix ("Confusion_matrix" function) and the mean average precision of the precision recall curve upon all possible queries.

In both sections it is possible to "uncomment" the ALLFEATPCA terms in order to project the descriptor "n"
space onto a lower dimentional space. In fact this trasformation is achieved by implementing the PCA 
(Eigen_Build.m) and then the "decriptor_projection.m" function, which reduces the descriptor space to the value 
of variable "dim". Note that you will have to call ALLFEATPCA, and not ALLFEAT, the terms inside the 
distance measure loop .
  

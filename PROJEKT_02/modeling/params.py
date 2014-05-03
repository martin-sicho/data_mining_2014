# NB CLASSIFICATION #
# crossvalidation params
NB_FOLDS = 10

# Support Vector Regression (SVR) #
VALIDATION_SET_SIZE = 0.0
SVR_FOLDS = 5
C = 5
EPSILON = 0.1
KERNEL = 'poly' # possible values: 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or a callable
                  # (more at http://scikit-learn.org/stable/modules/svm.html#svm-kernels)
DEGREE = 2 # only significant in poly, rbf, sigmoid
GAMMA = 0.0 # kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features will be taken
COEF0 = 1.5 # independent term in kernel function. It is only significant in poly/sigmoid
PROBABILITY = False # Whether to enable probability estimates. This must be enabled prior to calling fit, and will slow down that method.
VERBOSE = False

The files in this folder contain the prediction results for all sets included in the analysis. 

The filename format is as follows: 
[dataset_name]_prediction_results_[gt/at]_[model].csv

The gt/at is whether the model was trained on the whole of the training set (at) or only on the portion of the training set which has good templates (gt, see methods).

The model portion shows which feature sets were used for training and prediction. 

Results from 10-fold CV on the STCRDab set are also included.
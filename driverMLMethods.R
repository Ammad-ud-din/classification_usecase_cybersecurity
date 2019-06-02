rm(list=ls())
gc()
CURRENT_PATH = './classification_usecase_cybersecurity'
setwd(CURRENT_PATH)
SOURCE_PATH = './source_code/'
R_DATA_PATH = './Rdata/'
source(sprintf('%s/helpingfunctions.R',SOURCE_PATH))


#Run Cross Validation experiments on the training data
RESULT_PATH = './results/ml_methods_trainingdata/'
runMLMethodsTrainingData(R_DATA_PATH,RESULT_PATH)

collectMLMethodResultsTraininData(CURRENT_PATH,RESULT_PATH)

dimRed = "PCA"
sam = "down"
randomSeed = 5160
load(sprintf('%s/results/MLMethods_CV_Performances_%s_%s_Seed_%d.RData',CURRENT_PATH,dimRed,sam,randomSeed))
cvPerformances.mat
cvPerformances.cMat

sam = "orig"
load(sprintf('%s/results/MLMethods_CV_Performances_%s_%s_Seed_%d.RData',CURRENT_PATH,dimRed,sam,randomSeed))
cvPerformances.mat
cvPerformances.cMat

#Run ML Model on the holdout data
RESULT_PATH = './results/ml_methods_holdout_data/'
resfile = runMLMethodHoldoutData(R_DATA_PATH,RESULT_PATH)
load(resfile)
print(accuracy)

#Run ML Model on test data
RESULT_PATH = './results/test_set_predictions/'
resFile=runMLMethodTestData(R_DATA_PATH,RESULT_PATH)

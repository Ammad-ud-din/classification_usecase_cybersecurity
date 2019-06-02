loadLibraries <- function()
{
  library(ComplexHeatmap)
  library("circlize")
  library('gtools')
  library('gplots')
  library('RColorBrewer')
  library('R.matlab')
  library(reshape2)
  library(grid)
  library(gridExtra)
  library(reshape)
  library('Rmisc')
  library(ggthemes)
  library(ggplot2)
  library(grid)
  
}
theme_Publication <- function(base_size=14, base_family="Times") {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.text.x = element_text(angle = 45, hjust = 1,size = rel(1)),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size = unit(1, "cm"),
           legend.key.width = unit(2.5, "cm"),
           legend.margin = margin(t=1,b=1,r=1,l=1,unit = 'mm'),
           legend.title = element_blank(),#element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  library(RColorBrewer)
  col <- brewer.pal(n = 12, name = "Set3")
  discrete_scale("fill","Publication",manual_pal(values = c(col,"#000000")), ...)
  #c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","#000000")
}

scale_colour_Publication <- function(...){
  library(scales)
  library(RColorBrewer)
  col <- brewer.pal(n = 12, name = "Set3")
  #tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
  #c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C")
  
  discrete_scale("colour","Publication",manual_pal(values = c(col,"#000000")), ...)
  
}
or <- function(x){  as.numeric(sum(x!=0)>0) }

#####################################
###### Sorting Functions ############
#####################################
my.sort <- function(x,d) { sort(x,decreasing=d) }

#####################################
###### NA Functions #################
#####################################
mean.na <- function(x) { mean(na.omit(x)) }
sd.na <- function(x) { sd(na.omit(x)) }
min.na <- function(x) { min(na.omit(x)) }
max.na <- function(x) { max(na.omit(x)) }
sum.na <- function(x) { sum(na.omit(x)) }
which.na <- function(x) { which(na.omit(x)) }
length.na <- function(x) { sum(!is.na(x)) }
sort.na.dec <- function(x) { my.sort(na.omit(x),d=TRUE) }
sort.na.inc <- function(x) { my.sort(na.omit(x),d=FALSE) }
na.omit.print <- function(x) { a = na.omit(x); print(a[1:nrow(a),]) }

stratified_kfolds <- function(data,k){
  
  data = data[,order(apply(data,2,sum))]
  if(k>sum(data[,1])) warning("Number of folds should not be greater than smallest class size")
  
  inds = which(data[,1] == 1)
  finds = sample(length(inds),length(inds),replace=FALSE); #randomize the sequence 
  foldid = ceiling(finds/(length(inds)/k))
  names(foldid) = inds
  for(d in 2:ncol(data))
  {
    inds = which(data[,d] == 1)
    inds = inds[!(inds %in% names(foldid))] #remove already assigned samples
    finds = sample(length(inds),length(inds),replace=FALSE); #randomize the sequence 
    tmp = ceiling(finds/(length(inds)/k))
    names(tmp) = inds
    foldid = c(foldid,tmp)
  }
  return(foldid)
}
get_stratified_kfolds = function(Y=NULL,randomSeed=NULL,K=10)
{
  Ytrain.mat=matrix(0,nrow(Y),ncol=2)
  ind = which(Y[,1]==1)
  Ytrain.mat[ind,1]=1
  ind = which(Y[,1]==(-1))
  Ytrain.mat[ind,2]=1
  
  colSums(Ytrain.mat)
  colnames(Ytrain.mat)=c("1","-1")
  #K=10
  set.seed(randomSeed)
  foldid.10 = stratified_kfolds(data= Ytrain.mat,K)
  #table(foldid.10)
  return(list(foldid.10=foldid.10))
  
}

get_holdout_stratified_kfolds = function(Y=NULL,X=NULL,randomSeed=NULL,K=10)
{
  Ytrain.mat=matrix(0,nrow(Y),ncol=2)
  ind = which(Y[,1]==1)
  Ytrain.mat[ind,1]=1
  ind = which(Y[,1]==(-1))
  Ytrain.mat[ind,2]=1
  
  colSums(Ytrain.mat)
  colnames(Ytrain.mat)=c("1","-1")
  nonBinCols = colnames(Ytrain.mat)
  #K=10
  set.seed(randomSeed)
  foldid = stratified_kfolds(Ytrain.mat,K)
  #Y.holdout = Ytrain[which(foldid == 1),]
  Y.holdout = as.matrix(Y[which(foldid == 1),1])
  colSums(Y.holdout)
  X.holdout = X[which(foldid == 1),]
  
  Ytrain.mat = Ytrain.mat[-which(foldid == 1),]
  Xtrain = X[-which(foldid == 1),]
  Ytrain = as.matrix(Y[-which(foldid == 1),1])
  #foldid.2 = stratified_kfolds(data= Y.holdout[,!(colnames(Y.holdout)%in%nonBinCols)],2)
  set.seed(randomSeed)
  foldid.10 = stratified_kfolds(data= Ytrain.mat,K)
  #table(foldid.10)
  return(list(X=Xtrain,Y=Ytrain,X.holdout=X.holdout,Y.holdout=Y.holdout,foldid.10=foldid.10))

}

loadData <- function(DATA_PATH=NULL,trainingDataFile='train_data', trainingDataLabelFile='train_labels',testDataFile='test_data')
{
  trainingData = read.csv(sprintf("%s/%s.csv",DATA_PATH,trainingDataFile),header = F)
  #trainingData = apply(trainingData,2,null_to_na)
  trainingDataLabels = read.csv(sprintf("%s/%s.csv",DATA_PATH,trainingDataLabelFile),header = F)
  
  testData = read.csv(sprintf("%s/%s.csv",DATA_PATH,testDataFile),header = F)
  
  return(list(trainingData=trainingData,trainingDataLabels=trainingDataLabels,testData=testData))
}

getCustomRF = function()
{
  library(randomForest)
  customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
  customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
  customRF$grid <- function(x, y, len = NULL, search = "grid") {}
  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
  }
  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)
  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")
  customRF$sort <- function(x) x[order(x[,1]),]
  customRF$levels <- function(x) x$classes
  
  return(customRF)
}

runMLMethodsTrainingData = function(R_DATA_PATH=NULL,
                                    RESULT_PATH=NULL)
{

  #load required libraries/packages
  library(caret)
  library(glmnet)  
  
  
  #
  #source(sprintf('%s/helpingfunctions.R',SOURCE_PATH))
  #fix seed to generate random seeds, in case multiple cv runs are needed
  set.seed(1606)
  randomSeed = sample(100:10000,size=5,replace = F)[1]
  #number of cv folds
  numOfFolds = 10
  #load already saved RData obeject containing into training, test, holdout and cv foldids
  load(sprintf("%s/cv_train_labels_test_seed_%d.RData",R_DATA_PATH,randomSeed))
  
  #setting the number of samples
  numOfSamples = dim(X)[1]
  numOfFeatures = dim(X)[2]
  #dimensioality reduction with PCA
  dimRed = "PCA"
  #ML methods to be tested using CV
  #RF=Random Forest
  #DR=Decision Trees
  #KNN=K nearest neighbors
  #LR=Logistic linear regression
  MLMethods = c('RF','DT','KNN','LR')
  
  #Sampling approach to use
  # down = down sampling the majority class (to tackle class imbalances)
  # orig = no sampling is done
  samplingMethods = c("down","orig")
  
  for(sam in 1:length(samplingMethods))
  {
    
    for(m in 1:length(MLMethods))
    {
      #setting the file path where results will be saved
      resFile = sprintf('%s/CV_Performance_%s_%s_%s_%d.RData',RESULT_PATH,MLMethods[m],samplingMethods[sam],dimRed,randomSeed)
      
      #checking if the result already exists, then skip run this iteration
      if(file.exists(resFile))
      {
        print("Results Already Exists")
      }else{
        
        #Creating object to save predictions
        YPred = rep(NA,length(Y))
        #converting the Y to a character vector, in case it should be not be mislead with regression
        #since most model implementations have both regression and classification calls with a same function
        YObs = Y
        YObs[YObs==1]='A'
        YObs[YObs==(-1)]='B'
        #loop over CV folds
        for(fold in 1:numOfFolds){
          
          #getting training and test set indexes
          trainingIndex = which(fold!=foldid.10)
          testIndex = which(fold==foldid.10)
          
          #partition the data into training and test set, for the current CV-fold
          Xtrain.cv = as.matrix(X[trainingIndex,])
          Xtest.cv = as.matrix(X[testIndex,])
          
          Ytrain.cv = factor(YObs[trainingIndex], levels=unique(as.vector(YObs)))
          
          
          #getting the stratified fold ids for the nested CV set-up
          nestedCVData = get_stratified_kfolds(Y = as.matrix(Y[trainingIndex,1]),randomSeed,numOfFolds)
          
          folds.list.test = lapply(1:numOfFolds, function(x) as.vector(which(nestedCVData$foldid.10==x)))
          folds.list.train = lapply(1:numOfFolds, function(x) as.vector(which(nestedCVData$foldid.10!=x)))
          
          #adjusting the broader level experimental settings for the all models 
          if(sam==2)
          {
            #in case no sampling is done
            trctrl <- caret::trainControl(method = "cv", number = 10,
                                          savePredictions = TRUE,
                                          classProbs = TRUE,
                                          index = folds.list.train,
                                          indexOut = folds.list.test)
            
          }else{
            #in case "down" sampling is required
            trctrl <- caret::trainControl(method = "cv", number = 10,
                                          savePredictions = TRUE,
                                          classProbs = TRUE,
                                          sampling = "down",
                                          index = folds.list.train,
                                          indexOut = folds.list.test)
            
          }
          
          
          #setting the path where output of the PCA is stored. 
          #This is needed so that PCA should be run once for each CV fold 
          #it should not be repeated for each model. Saves computation time!
          foldPaths = sprintf('%s/Fold_%d_%s_%d_Var2.RData',R_DATA_PATH,fold,dimRed,randomSeed)
          
          
          if(file.exists(foldPaths))
          {
            #if PCA is already done for the current fold, load it directly
            load(foldPaths)
            print(dim(Xtrain.cv))
            
          }else{
            
            #if PCA hasn't been done, 
            
            # compute feature wise mean from the training data
            mFeatures = apply(Xtrain.cv,2,mean)
            
            # data centering: subtracting the mean from the training and test data
            Xtrain.cv = Xtrain.cv - matrix(as.vector(mFeatures),nrow=length(trainingIndex),ncol=ncol(Xtrain.cv))
            Xtest.cv = Xtest.cv - matrix(as.vector(mFeatures),nrow=length(testIndex),ncol=ncol(Xtest.cv))
            
            
            #Apply PCA on the training data
            #center = FALSE: since data is manually centered 
            #scale = FALSE: although scaling the data is also recommended but 
            #in this problem I would like that PCA components should capture interesting variation 
            #if I scale the data, then I will lose it   
            Xtrain.cv = prcomp(as.matrix(Xtrain.cv),center = FALSE,scale. = FALSE)
            Xtest.cv = predict(Xtrain.cv, newdata = as.matrix(Xtest.cv))
            
            #automatically select the components that explains 0.5% or more of the variation in the data.
            #this is also evident from the data exploration 
            pVar = sum(Xtrain.cv$sdev) * 0.005
            nPcs = which(Xtrain.cv$sdev >= pVar)
            Xtrain.cv = Xtrain.cv$x[,nPcs]
            Xtest.cv = Xtest.cv[,nPcs]
            print(dim(Xtrain.cv))
            #save the PCA object
            save(Xtrain.cv,Xtest.cv,file=foldPaths)
            
          }
          
          #run sparse linear logistic regression
          if(MLMethods[m]=='LR')
          {
            #set the parameters choices to be tuned with nested CV
            #param alpha: regularization alpha=0 (ridge), alpha=1 (lasso), alpha=[0.1..0.9] (elastic net)
            #param lambda: penality
            grid <- expand.grid(alpha = c(seq(0,1,by=0.1)), 
                                lambda =  exp(seq(log(0.001), 
                                                  log(5), length.out=100)))
            set.seed(1606)
            #perform the nested CV and train the model with CV parameter values
            model <- train(x=Xtrain.cv,y=Ytrain.cv, method = "glmnet",
                           trControl=trctrl,
                           tuneGrid = grid,
                           metric = 'ROC')
            #use the optimal model and parameters (that provided the best accuracy) to predict the test set.
            ypred.cv = predict(model, newdata = Xtest.cv)
          }
          if(MLMethods[m]=='KNN')
          {
            #KNN: parameter to be tuned  k = number of nearest neighbors
            set.seed(1606)
            model <- train(x=Xtrain.cv,y=Ytrain.cv, method = "knn",
                           trControl=trctrl,
                           tuneGrid = expand.grid(.k=1:30),
                           metric = 'ROC')
            
            
            ypred.cv = predict(model, newdata = Xtest.cv)
            
          }
          
          if(MLMethods[m]=='DT')
          {
            # let the function to tune parameters from the default settings
            set.seed(1606)
            model <- train(x=Xtrain.cv,y=Ytrain.cv, method = "rpart",
                           parms = list(split = "information"),
                           trControl=trctrl,
                           metric = 'ROC',
                           tuneLength = 10)
            ypred.cv = predict(model, newdata = Xtest.cv)
            
          }
          if(MLMethods[m]=='RF')
          {
            
            # defining the random forest function with custom settings
            # to set the parameter choices that need to be evaluated
            customRF = getCustomRF()
            
            
            tunegrid <- expand.grid(.mtry=c(1:ncol(Xtrain.cv)),
                                    .ntree=c(25,50,75,100,500,1000))
            
            
            
            set.seed(1606)
            model <- train(x=Xtrain.cv,y=Ytrain.cv, method = customRF,
                           trControl=trctrl,
                           tuneGrid = tunegrid,
                           metric = 'ROC')
            
            ypred.cv = predict(model, newdata = Xtest.cv)
            
          }
          
          
          #storing the current test set predictions
          YPred[testIndex] = as.character(ypred.cv)
          print(paste0('fold ',fold))
        }
        
        #once the CV round is completed, 
        #compute the accuracies including the confusion matrix and many more
        accuracy = confusionMatrix(as.factor(YPred),as.factor(YObs),positive='A')
        #put the predictions and observation into a list object, to be used later.
        predictions = list(YPred=YPred,YObs=YObs)
        
        #save the accuracies and prediction into a file. 
        save(accuracy,predictions,file=resFile)
        print(paste0("Saved: ",resFile))
      }
    }
  }

}

runMLMethodTestData = function(R_DATA_PATH=NULL,
                                    RESULT_PATH=NULL)
{
  library(caret)
  #see runMLMethodsTrainingData for inline comments
  set.seed(1606)
  randomSeed = sample(100:10000,size=5,replace = F)[1]
  numOfFolds = 10
  load(sprintf("%s/train_labels_test.RData",R_DATA_PATH,randomSeed))
  dimRed = "PCA"
  MLMethods = "RF"
  samplingMethods = "down"
  resFile = sprintf('%s/TestData_Predictions_%s_%s_%s_%d.RData',RESULT_PATH,MLMethods,samplingMethods,dimRed,randomSeed)
  if(file.exists(resFile))
  {
    print("Results Already Exists")
  }else{
    
    Y = data$trainingDataLabels
    YObs = Y
    YObs[YObs==1]='A'
    YObs[YObs==(-1)]='B'
    
    
    Xtrain = as.matrix(data$trainingData)
    Xtest = as.matrix(data$testData)
    
    Ytrain = factor(YObs[,1], levels=unique(as.vector(YObs[,1])))
    
    numOfSamples = dim(Xtrain)[1]
    numOfFeatures = dim(Xtest)[2]
    
    
    nestedCVData = get_stratified_kfolds(Y = as.matrix(Y),randomSeed,numOfFolds)
    
    folds.list.test = lapply(1:numOfFolds, function(x) as.vector(which(nestedCVData$foldid.10==x)))
    folds.list.train = lapply(1:numOfFolds, function(x) as.vector(which(nestedCVData$foldid.10!=x)))
    
    trctrl <- caret::trainControl(method = "cv", number = 10,
                                  savePredictions = TRUE,
                                  classProbs = TRUE,
                                  sampling = "down",
                                  index = folds.list.train,
                                  indexOut = folds.list.test)
    
    
    
    PCAOutputPath = sprintf('%s/TrainingData_%s_%d.RData',R_DATA_PATH,dimRed,randomSeed)
    
    if(file.exists(PCAOutputPath))
    {
      load(PCAOutputPath)
      print(dim(Xtrain))
      
      
    }else{
      
      
      mFeatures = apply(Xtrain,2,mean)
      
      Xtrain = Xtrain - matrix(as.vector(mFeatures),nrow=nrow(Xtrain),ncol=ncol(Xtrain))
      Xtest = Xtest - matrix(as.vector(mFeatures),nrow=nrow(Xtest),ncol=ncol(Xtest))
      
      Xtrain = prcomp(as.matrix(Xtrain),center = FALSE,scale. = FALSE)
      Xtest = predict(Xtrain, newdata = as.matrix(Xtest))
      
      pVar = sum(Xtrain$sdev) * 0.005
      nPcs = which(Xtrain$sdev >= pVar)
      Xtrain = Xtrain$x[,nPcs]
      Xtest = Xtest[,nPcs]
      print(dim(Xtrain))
      save(Xtrain,Xtest,file=PCAOutputPath)
      
    }
    
    
    
    customRF = getCustomRF()
    
    tunegrid <- expand.grid(.mtry=c(1:ncol(Xtrain)),
                            .ntree=c(25,50,75,100,500,1000))
    
    
    
    set.seed(1606)
    model <- train(x=Xtrain,y=Ytrain, method = customRF,
                   trControl=trctrl,
                   tuneGrid = tunegrid,
                   metric = 'ROC')
    
    ypred = predict(model, newdata = Xtest)
    
    predictions = as.character(ypred)
    predictions[predictions=='A']=1
    predictions[predictions=='B']=-1
    
    
    save(predictions,file=resFile)
    write.table(predictions,file=sprintf('%s/test_labels.csv',RESULT_PATH),sep = ",",row.names = FALSE,col.names = FALSE)
  }
  return(resFile)
}

runMLMethodHoldoutData = function(R_DATA_PATH=NULL,
                               RESULT_PATH=NULL)
{
  #see runMLMethodsTrainingData for inline comments
  library(caret)
  set.seed(1606)
  randomSeed = sample(100:10000,size=5,replace = F)[1]
  numOfFolds = 10
  load(sprintf("%s/cv_train_labels_test_seed_%d.RData",R_DATA_PATH,randomSeed))
  numOfSamples = dim(X)[1]
  numOfFeatures = dim(X)[2]
  dimRed = "PCA"
  MLMethods = "RF" 
  samplingMethods = "down"
  resFile = sprintf('%s/Holdout_Performance_%s_%s_%s_%d.RData',RESULT_PATH,MLMethods,samplingMethods,dimRed,randomSeed)
  
  if(file.exists(resFile))
  {
    print("Results Already Exists")
  }else{
    YObs.holdout = Y.holdout
    YObs.holdout[YObs.holdout==1]='A'
    YObs.holdout[YObs.holdout==(-1)]='B'
    
    YObs = Y
    YObs[YObs==1]='A'
    YObs[YObs==(-1)]='B'
    
    
    Xtrain.cv = as.matrix(X)
    Xtest.cv = as.matrix(X.holdout)
    
    Ytrain.cv = factor(YObs, levels=unique(as.vector(YObs)))
    
    nestedCVData = get_stratified_kfolds(Y = as.matrix(Y),randomSeed,numOfFolds)
    
    folds.list.test = lapply(1:numOfFolds, function(x) as.vector(which(nestedCVData$foldid.10==x)))
    folds.list.train = lapply(1:numOfFolds, function(x) as.vector(which(nestedCVData$foldid.10!=x)))
    
    
    trctrl <- caret::trainControl(method = "cv", number = 10,
                                  savePredictions = TRUE,
                                  classProbs = TRUE,
                                  sampling = samplingMethods,
                                  index = folds.list.train,
                                  indexOut = folds.list.test)
    
    
    PCAOutputPath = sprintf('%s/Holdout_%s_%d.RData',R_DATA_PATH,dimRed,randomSeed)
    
    
    if(file.exists(PCAOutputPath))
    {
      load(PCAOutputPath)
      print(dim(Xtrain.cv))
      
    }else{
      
      
      mFeatures = apply(Xtrain.cv,2,mean)
      
      Xtrain.cv = Xtrain.cv - matrix(as.vector(mFeatures),nrow=nrow(Xtrain.cv),ncol=ncol(Xtrain.cv))
      Xtest.cv = Xtest.cv - matrix(as.vector(mFeatures),nrow=nrow(Xtest.cv),ncol=ncol(Xtest.cv))
      
      Xtrain.cv = prcomp(as.matrix(Xtrain.cv),center = FALSE,scale. = FALSE)
      Xtest.cv = predict(Xtrain.cv, newdata = as.matrix(Xtest.cv))
      
      pVar = sum(Xtrain.cv$sdev) * 0.005
      nPcs = which(Xtrain.cv$sdev >= pVar)
      Xtrain.cv = Xtrain.cv$x[,nPcs]
      Xtest.cv = Xtest.cv[,nPcs]
      print(dim(Xtrain.cv))
      save(Xtrain.cv,Xtest.cv,file=PCAOutputPath)
      
    }
    customRF = getCustomRF()
    
    tunegrid <- expand.grid(.mtry=c(1:ncol(Xtrain.cv)),
                            .ntree=c(25,50,75,100,500,1000))
    
    
    
    set.seed(1606)
    model <- train(x=Xtrain.cv,y=Ytrain.cv, method = customRF,
                   trControl=trctrl,
                   tuneGrid = tunegrid,
                   metric = 'ROC')
    
    ypred = predict(model, newdata = Xtest.cv)
    
    
    
    
    YPred.holdout = as.character(ypred)
    
    accuracy = confusionMatrix(as.factor(YPred.holdout),as.factor(YObs.holdout),positive='A')
    predictions = list(YPred=YPred.holdout,YObs=YObs.holdout)
    
    save(accuracy,predictions,file=resFile)
    print(accuracy)
    print(paste0("Saved: ",resFile))
  }
  return(resFile)
}

collectMLMethodResultsTraininData = function(CURRENT_PATH=NULL,
                                  RESULT_PATH=NULL)
{
  set.seed(1606)
  randomSeed = sample(100:10000,size=5,replace = F)[1]
  dimRed = "PCA"
  MLMethods = c('RF','DT','KNN','LR')
  samplingMethods = c("down","orig")
  for(sam in samplingMethods)
  {
    cvPerformances.mat = matrix(NA,nrow=length(MLMethods),ncol = 5)
    rownames(cvPerformances.mat)=MLMethods
    colnames(cvPerformances.mat)=c('Balanced Accuracy','Sensitivity','Specificity','Precision','Recall')
    cvPerformances.cMat = vector(mode='list',length=length(MLMethods))
    names(cvPerformances.cMat)=MLMethods
    
    for(m in 1:length(MLMethods))
    {
      
      resFile = sprintf('%s/CV_Performance_%s_%s_%s_%d.RData',RESULT_PATH,MLMethods[m],sam,dimRed,randomSeed)
      
      if(file.exists(resFile))
      {
        load(resFile)
        cvPerformances.mat[m,]=as.numeric(accuracy$byClass[c(11,1,2,5,6)])
        cvPerformances.cMat[[m]]=accuracy$table
      }
    }
    save(cvPerformances.mat,cvPerformances.cMat,file=sprintf('%s/results/MLMethods_CV_Performances_%s_%s_Seed_%d.RData',CURRENT_PATH,dimRed,sam,randomSeed))
  }
}

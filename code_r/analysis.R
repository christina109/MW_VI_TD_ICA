# Author:  Christina Jin
# Contact: cyj.sci@gmail.com

rm(list = ls())

drive <- 'c:'

#f_sep <- '\\'
#drive  <- paste0(drive, ':', f_sep)


f_main <- file.path(drive, 'topic_mind wandering',  '3data2')
f_code <- file.path(f_main, 'code_r')
f_tool <- file.path(drive, 'files', 'tool_code')
f_result <- file.path(drive, 'topic_mind wandering', '4result2')

setwd(f_main)

source(file.path(f_code, 'feat_modelling2.r'))
source(file.path(f_code, 'analyze_beh.r'))
source(file.path(f_code, 'plot_data.r'))
source(file.path(f_tool, 'my func.r'))
library(beepr)
library(tcltk)
library(lme4)
library(BayesFactor)



###   global pars   ###
# subs
subs <- 301:330

#states.train <- c('nf', 'f')
states.train <- c()  # leave it  c() if use the same states set for both training and testing data sample
#states.train <- c('nf', 'f')
# states <- c('eoa', 'ioa')
#states       <- c('nf', 'f')  # must match with states.train; c(neg, pos)
states <- c('ot', 'mw')
probeLabOn.train <- TRUE
probeLabOn.test  <- TRUE
timeLabOn.train  <- FALSE
timeLabOn.test   <- FALSE
nBack  <- 3 
trainNbackOn <- TRUE
testNbackOn <- TRUE

# tasks
tasks  <- c('vs')        # tasks to build models
#test.tasks <- c('sart')  # tasks to test models, set to c() if no test tasks
test.tasks <- c('vs', 'sart')

# marker directory
f_measure_matfile <- 'measure_matfile'

# filter features
# leave default to load all feats
#filter.feats(list(c('', 'power'), c('', 'ISPC'), c('','ERP')))
#filter.feats()
filter.feats(list(), paste0('alphaEOICo', 1:32))


################################
###   MODELLING & PLOTTING   ###
################################

###   overall modelling   ###
# each data is a testing sample in each loop

# before: load raw data from all subs, combine data from both tasks, normalize

get.settings()

# set ML pars
mtype     <- 'svm'  # algorithm: lr, svm, knn
cPowerVec <- -5:15
gPowerVec <- -15:3
validType <- 'cv'

f_parMat  <- c()  # leave it to c() if no parMat will be loaded
gridSearchOn  <- FALSE

# normalize training / testing ?
if (any(unlist(lapply(measures, function(x){grepl('PC', x) || grepl('IC', x)})) && 
        !unlist(lapply(folders, function(x){grepl('epochs', x) || grepl('EO', x)})))) {
  print('RAW PCA features included. Turn OFF normalization.')
  normalizeOn <- FALSE
} else {
  print('Normalization is ON.')
  normalizeOn <- TRUE
}


# intialize
accMat <- matrix(0, length(subs), length(tasks) * (length(test.tasks) + 2))
colnames(accMat) <- c(paste0(tasks, '.cv'), paste0(tasks, '.loso'), paste0( tasks, '.test.', test.tasks))
kapMat <- accMat; senMat <- accMat; speMat <- accMat

n <- length(subs) * length(tasks)
i <- 0
pb <- tkProgressBar(title = 'LOPOCV', min = 0, max = n, width = 300)

for (taski in 1:length(tasks)){
  
  if (length(f_parMat) > 0) {
    parList = list(c = parMat[taski,'c'], gamma = parMat[taski, 'gamma'])
  } else {
    parList = list()
  }
  task <- tasks[taski]
  
  for (subi in 1:length(subs)){
    
    print(paste0('Iteration ', subi + length(subs)*(taski-1), ' out of ', n))
    
    subs2test  <- subs[subi]
    subs2train <- subs[-subi]
    print(paste0('Training data are: ', paste0(subs2train, collapse = ' ')))
    print(paste0('Testing data are: ',paste0(subs2test, collapse = ' ')))
      
    train      <- get.all.data(subs2train, task, measures, feats, folders, probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack = ifelse(trainNbackOn, nBack, NaN), states = states, states.train = states.train, normOn = normalizeOn)
    test.os    <- get.all.data( subs2test, task, measures, feats, folders, probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack = ifelse(trainNbackOn, nBack, NaN), states = states, states.train = states.train, normOn = normalizeOn)
    
    # balance: might be not necessary for group training because the nearly equivalent amount of observations of both classes
    train.balan <- balance.class(train)
    
    # validation
    print('Cross-validate...')
    if (validType == 'loo'){
      perf <-    leave.one.out(train, mtype, parList = parList)  # [default] balanced = 'copy' for the training split
    } else {
      perf <- cross.validation(train, mtype, parList = parList)  # [default] balanced = 'copy' for the training split
    }
    accMat[subi, paste0(task, '.cv')] <- perf$accuracy
    kapMat[subi, paste0(task, '.cv')] <- perf$kappa
    senMat[subi, paste0(task, '.cv')] <- perf$sensitivity
    speMat[subi, paste0(task, '.cv')] <- perf$specificity
    
    # test: leave-one-sub-out
    print('Test model on the leave-out data...')
    m <- feat.modeling(train.balan, test.os, mtype, parList = parList)
    perf <- measure.performance(m$predictions, m$observations)
    accMat[subi, paste0(task, '.loso')] <- perf$accuracy
    kapMat[subi, paste0(task, '.loso')] <- perf$kappa
    senMat[subi, paste0(task, '.loso')] <- perf$sensitivity
    speMat[subi, paste0(task, '.loso')] <- perf$specificity
    
    # test: self-reports
    for (testi in 1:length(test.tasks)) {
      
      test.task <- test.tasks[testi]
      test <- get.all.data(subs2test, test.task, measures, feats, folders, probeLabOn = probeLabOn.test, timeLabOn = timeLabOn.test, nBack = ifelse(testNbackOn, nBack, NaN), states = states, states.train = c(), normOn = normalizeOn)
      
      print(paste0('Test model on a different labeling of ', test.task))
      m <- feat.modeling(train.balan, test, mtype, parList = parList)
      perf <- measure.performance(m$predictions, m$observations)
      accMat[subi, paste0(task, '.test.', test.task)] <- perf$accuracy
      kapMat[subi, paste0(task, '.test.', test.task)] <- perf$kappa
      senMat[subi, paste0(task, '.test.', test.task)] <- perf$sensitivity
      speMat[subi, paste0(task, '.test.', test.task)] <- perf$specificity
      
    }  # loop over test.tasks
    
    # save the tempoerary matrixes
    print('Autosave...')
    save(file = 'temp_featModeling_pooledSub.rdata', accMat, kapMat, senMat, speMat)
    i <- i + 1
    print(paste0(round(i/n * 100, 1), "% done"))
    setTkProgressBar(pb, i, label = paste(round(i/n * 100, 1), '% done'))
    
  }  # loop over subs
}  # loop over tasks
close(pb)
  
acc.df         <- as.data.frame(accMat)
kappa.df       <- as.data.frame(kapMat)
sensitivity.df <- as.data.frame(senMat)
specificity.df <- as.data.frame(speMat)
acc.df$sub         <- subs
kappa.df$sub       <- subs
sensitivity.df$sub <- subs
specificity.df$sub <- subs

f_acc <- paste0('featModeling_poolSubs_', mtype, 
                '_', measures[1], '_', measures[length(measures)],
                '_by', ifelse(probeLabOn.train, 'Rating', ifelse(timeLabOn.train, 'Timecourse', 'EOAvsIOA')), 
                '_', min(subs), '_', max(subs), '_all3back.rdata')
save(file =  f_acc, acc.df, kappa.df, sensitivity.df, specificity.df)

summary(acc.df)
summary(kappa.df)
summary(sensitivity.df)
summary(specificity.df)



###   test features   ###
source(paste0(f_code, f_sep, 'feat_testing_main.r'))

get.settings()

# parallel processing setting
subs.bk <- subs
subsetid <- NaN  # Nan - no parallel processing
parallelProcId <- FALSE  # FALSE or a number
if (is.vector(subsetid) && !is.nan(subsetid)) {
  subs <- subs.bk[subsetid]
  print(paste0('Process subset: ', subsetid))
  print(paste0('Pallell processing id: ', parallelProcId))
}

# set pars
mtype     <- 'lr'  # algorithm: lr, svm, knn
validType <- 'cv'  # 'cv'
taskMerge     <- FALSE

# pars when mtype == 'svm'
useParMat    <- FALSE
f_parMat     <- FALSE  # FALSE or .rdata name
gridSearchOn <- FALSE
searchBy     <- ''  # 'acc', 'kappa', 'sen-spe'

# special treat on data before training
permutateOn   <- FALSE
groupOn       <- TRUE  # train on group data
lopocvOn      <- TRUE  # leave-one-participant-out cross-validation 

# feature combo testing
featComboOn   <- FALSE  # start from the best feature, add one by one until no accuracy increase

# normalize training / testing ?
if (any(unlist(lapply(measures, function(x){grepl('PC', x) || grepl('IC', x)})) && 
        !unlist(lapply(folders, function(x){grepl('epochs', x) || grepl('EO', x)})))) {
  print('RAW PCA features included. Turn OFF normalization.')
  normalizeOn <- FALSE
} else {
  print('Normalization is ON.')
  normalizeOn <- TRUE
}

if (mtype == 'svm' && !gridSearchOn && !is.logical(f_parMat)){
  load(f_parMat)
  if (nrow(parMat) == 30) { 
    parMat <- parMat[subs,] 
    sub.test <- acc.df[subs,]$sub
    }
  # check sub id
  if (!all.equal(sub.test,subs)) {warning('Check if the specified subs match the parMat!')}
} else {
  parMat <- NaN
}

if (is.logical(permutateOn) & !permutateOn){
  
  if (!featComboOn){
    f_temp <- featTesting(subs, tasks, mtype, validType, parMat, taskMerge, useParMat,
                          permutateOn, gridSearchOn, searchBy, parallelProcId, normalizeOn,
                          groupOn, lopocvOn)
    
    load(f_temp)
  }
 
} else if (is.numeric(permutateOn)){
  
  # intialize
  temp <- matrix(0, length(subs), permutateOn)
  for (mi in 1:length(measures)){
    eval(parse(text = paste0('acc.', measures[mi], ' <- temp')))
  }
  
  n <- permutateOn
  i <- 0
  pb0 <- tkProgressBar('Permutation testing ... ', min = 0, max = n, width = 300)
  
  for (pi in 1:permutateOn){
    featTesting(subs, tasks, mtype, validType, f_parMat, taskMerge, useParMat, permutateOn = TRUE)
    load('temp_featTesting.rdata')
    if (taskMerge){
      for (mi in 1:length(measures)){
        measure <- measures[mi]
        eval(parse(text = paste0('acc.', measure, '[,pi] <- accPerFeat[, measure]')))
      }
    } else {
      for (mi in 1:length(measures)){
        measure <- measures[mi]
        for (ti in 1:length(tasks)){
          if (ti == 1) {
            eval(parse(text = paste0('temp <- accPerFeat.', tasks[ti], '[,measure]')))
          } else {
            eval(parse(text = paste0('temp <- cbind(temp, accPerFeat.', tasks[ti], '[,measure])')))
          }
        }
        eval(parse(text = paste0('acc.', measure, '[,pi] <- rowMeans(temp)')))
      }
    }
    
    save(file = 'temp_permutation.rdata', list = paste0('acc.', measures))
    i <- i + 1
    setTkProgressBar(pb0, i, label = paste(round(i/n*100, 1), "% done")) 
  } 
  close(pb0)
}

f_feat <- paste0('featTesting', ifelse(length(tasks)>1 && taskMerge, '_taskMerge', ''), 
                 ifelse(groupOn,'_pooledSubs',''),
                 ifelse(lopocvOn, '_lopocv', ''),
                 '_', mtype,'_',measures[1], '_', measures[length(measures)], 
                 '_by', ifelse(probeLabOn.train, 'Rating', ifelse(timeLabOn.train, 'Timecourse', 'EOAvsIOA')),
                 '.rdata')
if (length(tasks) == 1||taskMerge){
  save(file = f_feat, accPerFeat)
} else {
  save(file = f_feat, list = paste0('accPerFeat.', tasks))
}

f_perm <- 'permutation_feats_EOAvsIOA.rdata'
save(file = f_perm, list = paste0('acc.', measures))



###   permutation test: feature testing   ###
f_perm <- 'permutation_feats_EOAvsIOA.rdata'
f_perf <- 'featTesting_taskMerge_12vs35.rdata'

load(f_perm)
load(f_perf)
percentage.df <- matrix(0, nrow(accPerFeat), length(measures))
colnames(percentage.df) <- measures
rownames(percentage.df) <- rownames(accPerFeat)
for (mi in 1:length(measures)){
  measure <- measures[mi]
  temp <- accPerFeat[, measure]
  eval(parse(text = paste0('temp <- cbind(temp, acc.', measure, ')')))
  for (si in 1:nrow(temp)){
    temp2 <- rank(temp[si,])/ncol(temp)
    percentage.df[si,mi] <- temp2[1]
  }
}

subs2exclude <- c(1,17)
percentage.df <- percentage.df[setdiff(rownames(percentage.df), subs2exclude),]

  
  
###   plot feature performance   ###
sortOn <- TRUE
sortBy <- 'mean' # 'mean', 'lowCI'
wholeModelOn <- FALSE
markSigOn <- TRUE

autoSaveOn <- FALSE

f_feat  = 'featTesting_pooledSubs_lopocv_svm_alphaEOICo1_alphaEOICo32_byEOAvsIOA.rdata'
#f_whole = 'featModelling_lr_12vs35_taskMerge.rdata'
load(f_feat)
#mu <- 0.5214
mu <- 0.5077

if (wholeModelOn) {
  load(f_whole)
  # subset whole modelling data
  if (nrow(acc.df) == 30) {
    acc.df <- acc.df[as.integer(rownames(accPerFeat)),]
  }
}

if (wholeModelOn) { 
  df <- data.frame(accPerFeat, whole = acc.df[, 'accMat'], sub = acc.df$sub)
  df.sum <- arrangeData2plot(df, se = 'ci')  # plot 95% ci
  df.sum$key <- factor(df.sum$key, levels = c(measures, 'whole'))
  df.sum <- arrange(df.sum, key)
  df.sum$name <- c(measureNames, 'whole model')
  df.sum$color<- c(measureColors, '#000000')
  df.sum$type <- c(measureTypes, 'whole model')
  df.sum$type <- factor(df.sum$type, labels = c('whole model', 'ERP', 'power', 'ISPC'), levels = c('whole model', 'ERP', 'power', 'ISPC'))
  shapes <- c(11, 15, 16, 17) # matching df.sum$type
} else {
  
  if (is.vector(accPerFeat)) {
    df <- data.frame(key = names(accPerFeat), val = accPerFeat)
  } else {
    df <- as.data.frame(accPerFeat)
    df <- na.omit(df)  # e.g., training idv models on rating
  }

  if (!grepl('pooled', f_feat) || grepl('lopocv', f_feat)) { # if NOT pooled data (individual modelling/pooled lococv)
    df$sub <- as.integer(rownames(df))
    df.sum <- arrangeData2plot(df, se = 'ci')
  } else { # if pooled data 
    df.sum <- df
    df.sum$mean <- df.sum$val
  }
  
  df.sum$key <- factor(df.sum$key, levels = measures)
  df.sum <- arrange(df.sum, key)
  df.sum$name <- measureNames
  df.sum$color <- measureColors
  df.sum$shape <- measureShapes
  df.sum$type <- measureTypes
  if (length(unique(measureTypes)) > 1) {
    df.sum$type <- factor(df.sum$type, labels = c('ERP', 'power', 'ISPC'), levels = c('ERP', 'power', 'ISPC'))
  }
  shapes <- unique(measureShapes)
}

if (sortOn){
  if (sortBy == 'mean') {
    sortedLabels <- arrange(df.sum, desc(mean))$name
  } else if (sortBy == 'lowCI') {
    sortedLabels <- arrange(df.sum, desc(mean - se))$name
  } else {
    print('Invalid sorting method! Use mean instead.')
    sortedLabels <- arrange(df.sum, desc(mean))$name
  }
  df.sum$name <- factor(df.sum$name, levels = sortedLabels, labels =  sortedLabels)
  df.sum <- arrange(df.sum, name)
} else {
  df.sum <- arrange(df.sum, type)
  df.sum$name <- factor(df.sum$name, levels = unique(df.sum$name), labels =  unique(df.sum$name))
}

# compare single feat perf to 0.5: modify df.sum
if (!grepl('pooled', f_feat) || grepl('lopocv', f_feat)) {

  df.sum$t <- 0 
  df.sum$p <- 0
  df.sum$sig <- ''
  for (feati in 1:length(feats)) {
    res <- t.test(df[, measures[feati]], mu = mu)
    df.sum[df.sum$key == measures[feati], c('t', 'p', 'sig')] <- c(round(res$statistic, 3), round(res$p.value, 3), ifelse(res$p.value < 0.001, '***', ifelse(res$p.value < 0.01, '**', ifelse(res$p.value < 0.05, '*', ''))))
  }
}

# generate ggplot
if (!grepl('pooled', f_feat) || grepl('lopocv', f_feat)) {
  p <- plot.point.flip(df.sum, 'mean', 'type', 'name', 'ci', ylim = c(0.45, 0.65), yintercep = mu, colorName = 'color', shape = shapes)
} else {
  p <- plot.point.flip(df.sum, 'mean', 'type', 'name', NaN, ylim = c(0.45, 0.65), yintercep = mu, colorName = 'color', shape = shapes)
}
p <- p + labs(y = 'Accuracy', x = 'EEG marker', color = 'Marker type', shape = 'Marker type')
if (markSigOn) {
  p <- p + scale_x_discrete(labels = paste0(df.sum$sig, df.sum$name))
}

if (autoSaveOn) {ggsave(paste0(sub('.rdata', '', f_feat),'_corrcl.png'), device = 'png', plot = p, path = f_result, dpi = 300, width = 20, height = 20, units = 'cm')}

# compare single feat perf to whole model perf
t2wm.res <- data.frame(marker = measureNames, t = 0, p = 0, sig = '', stringsAsFactors = FALSE)
for (feati in 1:length(feats)) {
  res <- t.test(df[, measures[feati]], df[,'whole'], paired = TRUE)
  t2wm.res[feati, 2:4] <- c(round(res$statistic, 3), round(res$p.value, 3), ifelse(res$p.value < 0.001, '***', ifelse(res$p.value < 0.01, '**', ifelse(res$p.value < 0.05, '*', ''))))
}
t2wm.res$marker <- factor(t2wm.res$marker, levels = measureNames, labels = measureNames)
t2wm.res <- arrange(t2wm.res, t)

# plot significantly predictive markers
get.settings()
f_feat
feats2plot <- get.sig.feats(df.sum)
plot.freq.change(subs, tasks, states, contentsets, feats2plot, normalizeOn = FALSE, 
                 f_output = sub('Testing', 'Sig', sub('.rdata', '', f_feat)))

# compare between different f_perfs
f_perf2 <- 'featTesting_taskMerge_searchByAcc_12vs35.rdata'
accPerFeat1 <- accPerFeat
load(f_perf2)
t2perf2.res <- data.frame(marker = measureNames, t = 0, p = 0, sig = '', stringsAsFactors = FALSE)
for (feati in 1:length(feats)) {
  res <- t.test(accPerFeat[, measures[feati]], accPerFeat1[, measures[feati]], paired = TRUE)
  t2perf2.res[feati, 2:4] <- c(round(res$statistic, 3), round(res$p.value, 3), ifelse(res$p.value < 0.001, '***', ifelse(res$p.value < 0.01, '**', ifelse(res$p.value < 0.05, '*', ''))))
}
t2perf2.res$marker <- factor(t2wm.res$marker, levels = measureNames, labels = measureNames)
t2perf2.res <- arrange(t2wm.res, t)


###   plot ML result: ACC   ###

library(ggpubr)
library(lsr)

datasumed = FALSE
#f_perf <- 'featModeling_poolSubs_svm_EOICo1_32_byEOAvsIOA_301_330.rdata'
f_perf <- 'featModeling_poolSubs_svm_alphaEOICo1_alphaEOICo32_byTimecourse_301_330_all3back.rdata'
ptitle  <- paste0('SVM classifier of ', 
                  ifelse(grepl('byRating', f_perf, ignore.case = TRUE), 'SELF-RATED', 'MANIPULATED'), 
                  ' attention')
load(f_perf)

if (grepl('subs', f_perf, ignore.case = TRUE)) {
  task2plot <- c('vs.cv', 'vs.loso', 'vs.test.vs', 'vs.test.sart')
  taskNames <- c('Cross-validation', 'Test on leave-out data', 
                 'Test on self-rated data: vs', 'Test on self-rated data: sart')
  tcolors <- c("#E6E6FA", "#984ea3", "#87CEFA", "#5F9EA0")
} else {
  task2plot <- c('vs', 'vs.test', 'sart.test')
  taskNames <- c('Cross-validation', 'Self-rated state: VS', 'Self-rated state: SART')
  tcolors <- c("#E6E6FA", "#984ea3", "#4682B4", "#386cb0")
}

tcolors <- tcolors[1:length(task2plot)]

# arrange data 
df <- arrangeData2plot(acc.df, sumData = datasumed)
df$sub <- factor(df$sub)
df$key <- factor(df$key, labels = task2plot, levels = task2plot)

if (!datasumed){
  p <- ggplot(df, aes(x = sub, y = val, fill = key)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.5) +
    scale_fill_manual(name = '', values = tcolors, labels = taskNames) +
    scale_x_discrete(labels = 1:length(subs)) +  # x-tick names
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = ptitle, x = 'Subject', y = 'Prediction accuracy') +
    geom_hline(aes(yintercept = 0.5), color = '#5F9EA0', size = 1, linetype = 2, alpha = 0.5, show.legend = FALSE) +
    get.mytheme()
} else {
  p <- plot.bar(df, 'mean',groupName = 'key', errorValName = 'se', ylim = c(0, 1), autoColor = FALSE) +
    scale_x_discrete(labels = toupper(taskNames)) +
    labs(title = '', x = '', y = 'Prediction accuracy') + 
    scale_fill_manual(name = '', values = tcolors, labels = taskNames)
  p <- ggplot(df, aes(x = key, y = mean, fill = key)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.5) +
    # scale_fill_manual(name = 'Training task-testing task', values = tcolors, labels = taskNames) +
    scale_x_discrete(labels = toupper(taskNames)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = '', x = 'Training task-testing task', y = 'Prediction accuracy') +
    get.mytheme()
}

# output the plot
ggsave(paste0('accBars_', sub('rdata', 'png', f_perf)), device = 'png', plot = p, path = f_result, dpi = 300, width = 40, height = 20, units = 'cm')

# to do t-test between acc and chance level:
chance <- 0.5214 
#chance <- 0.5077
task <- 'vs.test.sart'
t.test(acc.df[,task], mu = chance)  # t-test for each task
#cohensD(acc.df[,task], mu = chance)
acc.arrange <- arrangeData2plot(acc.df, sumData= FALSE)  
t.test(acc.arrange$val, mu = chance)  # t-test for all

# modelling perf comparison
f_perf2 <- 'featModelling_svm_parByAcc_EOAvsIOA.rdata'
acc.df1 <- acc.df
load(f_perf2)
# check
if (!all.equal(acc.df$sub, acc.df1$sub)) {
  print('Unmatched subs!')
}
task <- 'vs'
t.test(acc.df[,task], acc.df1[,task], paired = TRUE)  # t-test for each task
cohensD(acc.df[,task], acc.df1[,task], method = 'paired')


###   plot ML result: ACC, sensitivity, specificity   ###

library(ggpubr)

f_perf <- 'featModelling_svm_OTvsMW_byEOAvsIOA_301_329.rdata'

if (file_test('-f', f_perf)) {
  load(f_perf)
  errorOn <- FALSE
} else {
  errorOn <- TRUE
}

indications <- c('acc','sensitivity','specificity')

if (grepl('subs', f_perf, ignore.case = TRUE)) {
  task2plot <- c('vs.cv', 'vs.loso', 'vs.test.vs', 'vs.test.sart')
  taskNames <- c('Cross-validation', 'Test on leave-out data', 
                 'Test on self-rated data: vs', 'Test on self-rated data: sart')
  tcolors <- c("#E6E6FA", "#984ea3", "#87CEFA", "#5F9EA0")
} else {
  task2plot <- c('vs', 'vs.test', 'sart.test')
  taskNames <- c('Cross-validation', 'Self-rated state: VS', 'Self-rated state: SART')
  tcolors <- c("#E6E6FA", "#984ea3", "#4682B4", "#386cb0")
}

for (indi in 1:length(indications)){
  indication <- indications[indi]
  eval(parse(text = paste0('temp <-', indication, '.df')))

  temp <- arrangeData2plot(temp, sumData = FALSE)
  temp$indication <- indication
  if (!'sub' %in% colnames(temp)) {temp$sub <- subs2plot}
  if (indi == 1){
    df <- temp
  } else {
    df <- rbind(df, temp)
  }
}
df$sub <- factor(df$sub)

# create plot list
plist <- list()
for (taski in 1:length(task2plot)) {
  task <- task2plot[taski]
  dfSub <- subset(df, key == task)
  p <- plot.point(dfSub, xName = 'sub', yName = 'val', condName = 'indication', ylim = c(0,1))
  p <- p + scale_x_discrete(labels = 1:length(subs)) +  # x-axis tick names
    scale_shape_manual(name = 'Performance Measure', values = 15:17, labels = c('Accuracy', 'Sensitivity', 'Specificity')) +
    scale_color_manual(name = 'Performance Measure', values = c('#4B0082', colors[2], colors[1]), labels = c('Accuracy', 'Sensitivity', 'Specificity'))
  if (length(task2plot) > 1){
    if (taski == 1){
      p <- p + labs(title = taskNames[taski], x = 'Subject', y = '')
    } else {
      p <- p + labs(title = taskNames[taski], x = '', y = '') 
    }
  }
  eval(parse(text = paste0('plist[[', taski, ']] <- p')))
}

# show plots
if (length(plist) > 1){
  ggarrange(plotlist = plist, common.legend = TRUE, nrow = length(task2plot), ncol = 1)
} else {
  plist
}

ggsave(paste0('perfPoints_', sub('rdata', 'png', f_perf)), device = 'png', path = f_result, dpi = 300, width = 25, height = 30, units = 'cm')
if (errorOn) {stop('File doesn\'t exist!' )}


###   analyze pred acc sum   ###
mu <- 0.5
f_acc <- 'featModeling_poolSubs_svm_alphaEOICo2_alphaEOICo2_byEOAvsIOA_301_330.rdata'
load(f_acc)
for (i in 1:4) {
  print(t.test(acc.df[,i], mu = mu))
}

library(equivalence)
f_acc <- 'featModeling_poolSubs_svm_alphaEOICo2_alphaEOICo2_byRating_301_330.rdata'
load(f_acc)
acc.df0 <- acc.df
f2_acc <- 'featModeling_poolSubs_svm_alphaEOICo2_alphaEOICo2_byTimecourse_301_330.rdata'
load(f2_acc)
for (i in 1:4) {
  print(t.test(acc.df0[,i], acc.df[,i], paired = TRUE))
  print(tost(acc.df0[,i], acc.df[,i], paired = TRUE))
}

###   analyze feats tesing   ###
#task.demands <- c(10, 11, 19, 30, 28, 16, 23, 8, 18, 3, 4, 2, 15, 6, 7, 5, 9)
#vigilance <- c(12, 3, 2, 7, 1, 15, 26)
#mental.state <- c(24, 9, 21, 17, 20, 2, 4, 5)
task.demands <- c(3,4,2,15,6,7,5,9)
vigilance <- c(1,26)
mental.state <- c(17,2,4)
tv <- intersect(task.demands, vigilance)
tm <- intersect(task.demands, mental.state)
vm <- intersect(vigilance, mental.state)
tvm <- intersect(tv, tm)
tv2 <- setdiff(tv, tvm)
vm2 <- setdiff(vm, tvm)
tm2 <- setdiff(tm, tvm)
t <- setdiff(setdiff(task.demands,tv),tm)
v <- setdiff(setdiff(vigilance,tv),vm)
m <- setdiff(setdiff(mental.state, tm), vm)
all <- union(union(task.demands, vigilance),mental.state)
all <- sort(all)



###  output beh data with labs   ###
library(R.matlab)
for (subi in 1:length(subs)) {
  sub <- subs[subi]
  task <- 'vs'
  beh <- get.beh.data(sub, task)
  beh <- subset(beh, select = -c(rating.response, rating.rt))
  beh <- na.omit(beh)
  
  writeMat(file.path(f_main, 'beh_matfile', paste0(sub, '.mat')), 
           cond = beh$block, blockId = beh$blocks.thisN,
           corr = beh$resp.corr, rt = beh$resp.rt, dis2pr = beh$dis2pr,
           pr_resp = beh$pr.resp, urevent = beh$urevent, 
           fatigue = beh$fatigue, state = beh$state)
}



###   analyze ratings   ###
conditions <- c('Counting', 'Non-counting', 'SART')

# intialize
df.mean <- matrix(0, length(subs), length(conditions)+1)
colnames(df.mean) <- c('sub', conditions)
df.mean[,1] <- subs

for (si in 1:length(subs)) {
  sub <- subs[si]
  tp <- get.beh.data(sub, 'vs', 1)
  temp <- subset(tp, block == 'eoa' & pr.resp != 0)
  df.mean[si, 'Counting'] <- mean(temp$pr.resp)
  temp <- subset(tp, block == 'ioa' & pr.resp != 0)
  df.mean[si, 'Non-counting'] <- mean(temp$pr.resp)
  tp <- get.beh.data(sub, 'sart', 1)
  temp <- subset(tp, pr.resp != 0)
  df.mean[si, 'SART'] <- mean(temp$pr.resp)
} 
df.tasks <- df.mean

# plot
df <- arrangeData2plot(df.mean)
df$key <- factor(df$key, levels = c('Counting', 'Non.counting', 'SART'), labels = conditions)
p <- plot.bar(df, 'mean', NaN, 'key', 'se')
p <- p + labs(x = '', y = 'Mean rating')
p1 <- p 

conditions <- c('First-half VS', 'Second-half VS')

# intialize
df.mean <- matrix(0, length(subs), length(conditions)+1)
colnames(df.mean) <- c('sub', conditions)
df.mean[,1] <- subs

for (si in 1:length(subs)) {
  sub <- subs[si]
  tp <- get.beh.data(sub, 'vs', 1)
  temp <- subset(tp, fatigue == 'nf' & pr.resp != 0)
  df.mean[si, 'First-half VS'] <- mean(temp$pr.resp)
  temp <- subset(tp, fatigue == 'f' & pr.resp != 0)
  df.mean[si, 'Second-half VS'] <- mean(temp$pr.resp)
} 
df.timecourse <- df.mean

# plot
df <- arrangeData2plot(df.mean)
df$key <- factor(df$key, levels = c('First.half.VS', 'Second.half.VS'), labels = conditions)
p <- plot.bar(df, 'mean', NaN, 'key', 'se')
p <- p + labs(x = '', y = '')
p2 <- p

# save
g <- ggarrange(p1, p2, nrow = 1, ncol = 2)
ggsave('rating.tif', device = 'tiff', plot = g, path = f_result, dpi = 300, width = 30, height = 14, units = 'cm')

# statistics 
library(ez)
df <- arrangeData2plot(df.tasks, FALSE)
df$sub <- factor(df$sub, levels = subs)
df$key <- factor(df$key, levels = c('Counting', 'Non.counting', 'SART'), labels = c('Counting', 'Non-counting', 'SART'))
m <- ezANOVA(df, dv = val, wid = sub, within = key)
m$ANOVA
pairwise.t.test(df$val, df$key, paired = T, p.adjust.method="bonferroni")

t.test(df.timecourse[,2], df.timecourse[,3], paired = TRUE)
library(lsr)
cohensD(df.timecourse[,2], df.timecourse[,3], method = 'paired')



###   analyze rating removal rate   ###
# removal = rating of zero

conditions <- c('Nonzero', 'Zero')

# intialize
df.count <- matrix(0, length(subs), length(conditions)+2)
colnames(df.count) <- c('sub', conditions, 'Sum')
df.count[,1] <- subs

for (si in 1:length(subs)) {
  sub <- subs[si]
  for (task in c('vs', 'sart')) {
    tp <- get.beh.data(sub, task, 1)
    temp <- subset(tp, pr.resp != 0)
    df.count[si, 'Nonzero'] <- df.count[si, 'Nonzero'] + nrow(temp)
    temp <- subset(tp, pr.resp == 0)
    df.count[si, 'Zero'] <- df.count[si, 'Zero'] + nrow(temp)
  }
  df.count[si, 'Sum'] <- df.count[si, 'Nonzero'] + df.count[si, 'Zero']
}
  
removal <- df.count[,'Zero'] / df.count[,'Sum']
tp_sum <- colSums(df.count[,3:4])
tp_sum[1]/tp_sum[2]



###   analyze ACC & RTs   ###

logRTon <- TRUE

conditions <- c('Counting', 'Non-counting')

# intialize
df.rt <- matrix(0, length(subs), length(conditions)+1)
colnames(df.rt) <- c('sub', conditions)
df.rt[,1] <- subs
df.acc <- df.rt

for (si in 1:length(subs)) {
  sub <- subs[si]
  tp <- get.beh.data(sub, 'vs')
  temp <- subset(tp, block == 'eoa')
  df.acc[si, 'Counting'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(temp, resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
    }
  df.rt[si, 'Counting'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
  
  temp <- subset(tp, block == 'ioa')
  df.acc[si, 'Non-counting'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(temp, resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'Non-counting'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
} 
rt.demands <- df.rt
acc.demands <- df.acc

# plot
df <- arrangeData2plot(df.rt)
df$key <- factor(df$key, levels = c('Counting', 'Non.counting'), labels = conditions)
p1 <- plot.bar(df, 'mean', NaN, 'key', 'se')
df2 <- arrangeData2plot(df.acc)
df2$key <- factor(df2$key, levels = c('Counting', 'Non.counting'), labels = conditions)
p1 <- p1 + geom_point(data = df2, aes(x = key, y = mean*ifelse(logRTon,4,1)), 
                      shape = 1, size = 5, alpha = 1, stroke = 1.5, color = 'blue')
p1 <- p1 + geom_errorbar(data = df2, aes(ymin = (mean - se)*ifelse(logRTon,4,1), ymax = (mean + se)*ifelse(logRTon,4,1)),
                         width = 0.2, size = 0.8, color = 'blue')
p1 <- p1 + scale_y_continuous(limits = c(0, ifelse(logRTon,4,1)), sec.axis = sec_axis(~./ifelse(logRTon,4,1), name = 'Accuracy'))
p1 <- p1 + labs(x = '', y = ifelse(logRTon, 'Log-transformed response time (ms)', 'Mean response time (s)'), title = 'Visual Search')
p1 <- p1 + theme(axis.line.y.right = element_line(color = 'blue'), 
                 axis.ticks.y.right = element_line(color = 'blue'),
                 axis.text.y.right = element_text(color = 'blue'),
                 axis.title.y.right = element_text(color = 'blue'))

conditions <- c('First-half', 'Second-half')

# intialize
df.rt <- matrix(0, length(subs), length(conditions)+1)
colnames(df.rt) <- c('sub', conditions)
df.rt[,1] <- subs
df.acc <- df.rt

for (si in 1:length(subs)) {
  sub <- subs[si]
  tp <- get.beh.data(sub, 'vs')
  temp <- subset(tp, fatigue == 'nf')
  df.acc[si, 'First-half'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(temp, resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'First-half'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
  
  temp <- subset(tp, fatigue == 'f')
  df.acc[si, 'Second-half'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(temp, resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'Second-half'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
} 
rt.timecourse <- df.rt
acc.timecourse <- df.acc

# plot
df <- arrangeData2plot(df.rt)
df$key <- factor(df$key, levels = c('First.half', 'Second.half'), labels = conditions)
p2 <- plot.bar(df, 'mean', NaN, 'key', 'se')
df3 <- arrangeData2plot(df.acc)
df3$key <- factor(df3$key, levels = c('First.half', 'Second.half'), labels = conditions)
p2 <- p2 + geom_point(data = df3, aes(x = key, y = mean*ifelse(logRTon,4,1)), shape = 1, size = 5, alpha = 1, stroke = 1.5, color = 'blue')
p2 <- p2 + geom_errorbar(data = df3, aes(ymin = (mean -se)*ifelse(logRTon,4,1), ymax = (mean+se)*ifelse(logRTon,4,1)), width = 0.2, size = 0.8, color = 'blue')
p2 <- p2 + scale_y_continuous(limits = c(0, ifelse(logRTon,4,1)), sec.axis = sec_axis(~./ifelse(logRTon,4,1)))
p2 <- p2 + labs(x = '', y = '', title = 'Visual Search')
p2 <- p2 + theme(axis.line.y.right = element_line(color = 'blue'), 
                 axis.ticks.y.right = element_line(color = 'blue'),
                 axis.text.y.right = element_text(color = 'blue'),
                 axis.title.y.right = element_text(color = 'blue'))

conditions <- c('On-task', 'Mind-wandering')

# intialize
df.rt <- matrix(0, length(subs), length(conditions)+1)
colnames(df.rt) <- c('sub', conditions)
df.rt[,1] <- subs
df.acc <- df.rt
df.count <- df.acc

for (si in 1:length(subs)) {
  sub <- subs[si]
  tp <- get.beh.data(sub, 'vs', nBack)
  temp <- subset(tp, state == 'ot')
  df.count[si, 'On-task'] <- nrow(temp)
  df.acc[si, 'On-task'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(temp, resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'On-task'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
  
  temp <- subset(tp, state == 'mw')
  df.count[si, 'Mind-wandering'] <- nrow(temp)
  df.acc[si, 'Mind-wandering'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(temp, resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'Mind-wandering'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
} 
rt.reportsvs <- df.rt
acc.reportsvs <- df.acc
count.reportsvs <- df.count

# plot
df.rt <- na.omit(df.rt)
df <- arrangeData2plot(df.rt)
df$key <- factor(df$key, levels = c('On.task', 'Mind.wandering'), labels = conditions)
p3 <- plot.bar(df, 'mean', NaN, 'key', 'se')
df.acc <- na.omit(df.acc)
df4 <- arrangeData2plot(df.acc)
df4$key <- factor(df4$key, levels = c('On.task', 'Mind.wandering'), labels = conditions)
p3 <- p3 + geom_point(data = df4, aes(x = key, y = mean*ifelse(logRTon,4,1)), shape = 1, size = 5, alpha = 1, stroke = 1.5, color = 'blue')
p3 <- p3 + geom_errorbar(data = df4, aes(ymin = (mean + se)*ifelse(logRTon,4,1), ymax = (mean - se)*ifelse(logRTon,4,1)), width = 0.2, size = 0.8, color = 'blue')
p3 <- p3 + scale_y_continuous(limits=c(0,ifelse(logRTon,4,1)),sec.axis = sec_axis(~./ifelse(logRTon,4,1), name = ''))
p3 <- p3 + labs(x = '', y = '', title = 'Visual Search')
p3 <- p3 + theme(axis.line.y.right = element_line(color = 'blue'), 
                 axis.ticks.y.right = element_line(color = 'blue'),
                 axis.text.y.right = element_text(color = 'blue'),
                 axis.title.y.right = element_text(color = 'blue'))

conditions <- c('On-task', 'Mind-wandering')

# intialize
df.rt <- matrix(0, length(subs), length(conditions)+1)
colnames(df.rt) <- c('sub', conditions)
df.rt[,1] <- subs
df.acc <- df.rt
df.count <- df.acc

for (si in 1:length(subs)) {
  sub <- subs[si]
  tp <- get.beh.data(sub, 'sart', nBack)
  temp <- subset(tp, state == 'ot') # all accuracy
  df.count[si, 'On-task'] <- nrow(temp)
  #temp <- subset(tp, state == 'ot' & type == 'target')  # only "NOGO" accuracy
  df.acc[si, 'On-task'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(tp, state == 'ot' & resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'On-task'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
  
  temp <- subset(tp, state == 'mw') # all accuracy
  df.count[si, 'Mind-wandering'] <- nrow(temp)
  #temp <- subset(tp, state == 'mw' & type == 'target')  # only "NOGO" accuracy
  df.acc[si, 'Mind-wandering'] <- mean(temp$resp.corr, na.rm = TRUE)
  temp <- subset(tp, state == 'mw' & resp.corr == 1)
  rtarr <- temp$resp.rt
  if (logRTon) {
    rtarr <- rtarr * 1000 # turn s to ms 
    rtarr <- log10(rtarr)
  }
  df.rt[si, 'Mind-wandering'] <- mean(rtarr, na.rm = TRUE, trim = 0.025)
} 
rt.reportssart <- df.rt
acc.reportssart <- df.acc
count.reportssart <- df.count

# plot
df.rt <- na.omit(df.rt)
df <- arrangeData2plot(df.rt)
df$key <- factor(df$key, levels = c('On.task', 'Mind.wandering'), labels = conditions)
p4 <- plot.bar(df, 'mean', NaN, 'key', 'se')
df.acc <- na.omit(df.acc)
df5 <- arrangeData2plot(df.acc)
df5$key <- factor(df5$key, levels = c('On.task', 'Mind.wandering'), labels = conditions)
p4 <- p4 + geom_point(data = df5, aes(x = key, y = mean*ifelse(logRTon,4,1)), shape = 1, size = 5, alpha = 1, stroke = 1.5, color = 'blue')
p4 <- p4 + geom_errorbar(data = df5, aes(ymin = (mean-se)*ifelse(logRTon,4,1), ymax = (mean+se)*ifelse(logRTon,4,1)), width = 0.2, size = 0.8, color = 'blue')
p4 <- p4 + scale_y_continuous(limits= c(0, ifelse(logRTon,4,1)), sec.axis = sec_axis(~./ifelse(logRTon,4,1), name = ''))
p4 <- p4 + labs(x = '', y = '', title = 'SART')
p4 <- p4 + theme(axis.line.y.right = element_line(color = 'blue'), 
                 axis.ticks.y.right = element_line(color = 'blue'),
                 axis.text.y.right = element_text(color = 'blue'),
                 axis.title.y.right = element_text(color = 'blue'))

# save
g <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggsave(paste0('beh scores',ifelse(logRTon, '_log',''), '.tif'), device = 'tiff', plot = g, path = f_result, dpi = 300, width = 26, height = 24, units = 'cm')

# statistics 
library(lsr)
tp <- rt.reportsvs  # change this
tp <- na.omit(tp)
ttestBF(tp[,2]-tp[,3], paried = TRUE)


t.test(tp[,2], tp[,3], paired = TRUE)
cohensD(tp[,2], tp[,3], method = 'paired')

# nonnormality check 
shapiro.test(tp[,2])
shapiro.test(tp[,3])


# selection analysis
vsIn <- count.reportsvs[,'On-task'] >0 & count.reportsvs[,'Mind-wandering']
sartIn <- count.reportssart[,'On-task'] >0 & count.reportssart[,'Mind-wandering']
selection <- 1:30
selection <- paste0(ifelse(vsIn, '*', ''), selection)
selection <- paste0(ifelse(sartIn, '+', ''), selection)
df2write <- data.frame(selection, count.reportsvs[,2:3], count.reportssart[,2:3])



###   count prac trials   ###
tasks2count <- c('sart', 'vs')  # avoid tasks set in the global  
pracCount  <- matrix(0, length(subs), length(tasks2count))
colnames(pracCount) <- tasks2count
rownames(pracCount) <- subs
for (si in 1:length(subs)) {
  sub <- subs[si]
  for (ti in 1:length(tasks2count)) {
    task <- tasks2count[ti]
    temp <- get.prac.beh(sub, task)
    pracCount[si, ti] <- nrow(temp)
  }
}
colMeans(pracCount)
colSds(pracCount)


###   count dataset size   ###
tasks <- c('vs')
statesets <- list(c('ot','mw'), c('nf', 'f'), c('eoa', 'ioa'))
classifiers <- c('Rating', 'Timecourse', 'Task')
probeLabOn <- c(TRUE, FALSE,FALSE)
timeLabOn <- c(FALSE, TRUE, FALSE)
df <- data.frame(task = rep(tasks, each = length(classifiers)), 
                 classifier = classifiers, N = 0)
for (ti in 1:length(tasks)) {
  task <- tasks[ti]
  for (si in 1:3){
    states <- statesets[[si]]
    tp <- get.all.data(301:330, task, measures[1], feats[1], folders[1], 
                       probeLabOn = probeLabOn[si], timeLabOn = timeLabOn[si], 
                       nBack = ifelse(si==1, 3, NaN), 
                       states = states, states.train = c(), 
                       normOn = FALSE)
    df[(ti-1)*2+si, 1] = task
    df[(ti-1)*2+si, 2] = classifiers[si]
    df[(ti-1)*2+si, 3] = nrow(tp)
  }
}
save(file = 'datasize.rdata', df) # sts are from matlab 
library(R.matlab)
library(matrixStats)
library(e1071)
library(lme4)
library(class)
library(caret)
library(DMwR)
library(utils)
library(tcltk)
library(neuralnet)



get.settings <- function(){
  print(paste0('Current feature path: ', f_measure_matfile))
  print(paste0('Tasks to build models: ', paste0(tasks, collapse = ' & ')))
  print(paste0('Current states for training: ', paste0(states.train, collapse = ' vs. ')))
  print(paste0('labeled by ', ifelse(probeLabOn.train, 'self-reports', ifelse(timeLabOn.train, 'time course', 'exp manipulations'))))
  print('--------------------------------------------------')
  print(paste0('Tasks to test models: ', paste0(test.tasks, collapse = ' & ')))
  print(paste0('Current states of interest: ', paste0(states, collapse = ' vs. ')))
  print(paste0('labeled by ', ifelse(probeLabOn.test, 'self-reports', ifelse(timeLabOn.test, 'time course', 'exp manipulations'))))
  print('--------------------------------------------------')
  print(paste0('Subject count: ', length(subs)))
  print(paste0('Current EEG markers are ', paste0(measures, collapse = ', ')))
  
}



grid.search.svm <- function(data, cVec, gVec, validType = 'cv', balanced = 'copy', searchBy = 'acc', checkOn = FALSE) {
  
  mat <- matrix(0,length(cVec), length(gVec))
  count <- length(cVec) * length(gVec)
  ti <- 1
  pb <- tkProgressBar(title = 'Grid search SVM', min = 0, max = count, width = 500)
  for (ci in 1:length(cVec)){
    for (gi in 1:length(gVec)){
      parList <- list(c = cVec[ci], gamma = gVec[gi])
      
      # only for check purpose
      if (checkOn) {tic(paste0('Grid search when c=', parList$c, ', gamma=', parList$gamma))}
      
      if(validType == 'loo'){
        perf <- leave.one.out(data, mtype = 'svm', parList = parList, balanced = balanced)
      } else {
        perf <- cross.validation(data, mtype = 'svm', parList = parList, balanced = balanced)
      }
      
      # fail to build a model?
      if (is.nan(perf$accuracy)){
        print('No model being built. Quit grid search of the input data')
        close(pb)
        return(NaN)
      }
      
      if (searchBy == 'acc'){
        mat[ci,gi] <- perf$accuracy
      } else if (searchBy == 'kappa'){
        mat[ci,gi] <- perf$kappa
      } else if (searchBy %in% c('sen-spe', 'spe-sen')) {
        mat[ci,gi] <- perf$sensitivity + perf$specificity
      } else {
        warning('Invalid search criteria. Use the highest accuracy.')
        mat[ci,gi] <- perf$accuracy
      }
      
      if (checkOn){toc()}
      
      setTkProgressBar(pb, ti, label = paste(round(ti/count * 100, 1), "% done"))
      ti <- ti + 1
    }
  }
  close(pb)
  return(mat)
  
}



leave.one.out <- function(data, mtype, parList = list(), balanced = 'copy'){
  
  print('Model fitting. Use LEAVE-ONE-OUT cross-validation')
  for (i in 1:nrow(data)) {
    
    test  <- data[i, ]
    train <- data[-i,]
    
    if (is.character(balanced)) {
      train <- balance.class(train, balanced)
    }
    
    pred <- feat.modeling(train, test, mtype, parList)
    if (i == 1){
      pred_class <- pred$predictions
      obs_class  <- pred$observations
    } else {
      pred_class <- unlist(list(pred_class, pred$predictions))
      obs_class  <- unlist(list(obs_class, pred$observations))
    }
    
  }
  
  performance <- measure.performance(pred_class, obs_class)
  return(performance)
  
}



cross.validation <- function(data, mtype, nfold = 10, parList = list(), balanced = 'copy'){
  
  print(paste0('Model fitting. Use ', nfold, '-fold CROSS-VALIDATION'))
  
  # balanced indicates the way to balance the training dataset 
  
  idx <- split.sample(data, nfold)
  
  # unable split? (because insufficient samples)
  if (length(idx) == 0) {
    return(list(accuracy = NaN, kappa = NaN, 
                sensitivity = NaN, specificity = NaN))
  }
  
  for (foldi in 1:nfold) {
    
    train <- data[idx$trainIdx[[foldi]],]
    test  <- data[idx$testIdx[[foldi]], ]
    
    if (is.character(balanced)) {
      train <- balance.class(train, balanced)
    }
    
    pred <- feat.modeling(train, test, mtype, parList)
    if (foldi == 1){
      pred_class <- pred$predictions
      obs_class  <- pred$observations
    } else {
      pred_class <- unlist(list(pred_class, pred$predictions))  # to remain the factor type
      obs_class  <- unlist(list(obs_class, pred$observations))
    }
    
  }
  
  performance <- measure.performance(pred_class, obs_class)
  return(performance)
  
}



feat.modeling <- function(train, test, mtype, parList = list()){
  
  if (mtype == 'lr'){
    
    m <- glm(state ~ ., family = binomial('logit'), data = train)
    p <- predict(m, test, type = "response")
    states <- names(summary(train$state))
    if (is.factor(train$state)) {
      states <- factor(states, levels = states, labels = states)
    }
    p_class <- ifelse(p>0.5, states[2], states[1])
    p_class <- factor(p_class, levels = 1:length(states), labels = states)
    #pStrength <- p
    
  } else if (mtype == 'svm') {
    
    if (length(parList) > 0){
      m <- svm(state ~ ., train, probability = TRUE, cost = parList$c, gamma = parList$gamma)  
    } else {
      m <- svm(state ~ ., train, probability = TRUE)
    }

    p_class <- predict(m, test)
    #prob <- attr(predict(m, test, probability = TRUE), 'probabilities')
    #pStrength <- prob[,2] - prob[,1]
    
  } else if (mtype == 'knn') {
    
    # split feats from lab
    trainFeats <- train[, -which(colnames(train) %in% 'state')]
    testFeats  <-  test[, -which(colnames(test)  %in% 'state')]
    p_class <- knn(trainFeats, testFeats, train$state, k = 5)
    #p <- attr(knn(trainFeats, testFeats, train$state, k = 5, prob = TRUE),'prob')
    #p[p_class == states[1]] <- -p[p_class == states[1]]
    #pStrength <- p
    m <- 'No model for KNN'
    
  } else if (mtype == 'tree') {
    
    m <- rpart(state ~ ., train, method = "class", control = rpart.control(cp=0.001))
    p_class <- predict(m,test,type="class")
    
  } else if (mtype == 'nn') {
    temp <- names(train)
    train$state <- as.numeric(train$state) - 1
    f <- as.formula(paste("state ~", paste(temp[!temp %in% "state"], collapse = " + ")))
    nfeats <- length(temp[!temp %in% "state"])
    m <- neuralnet(f, data = train, hidden = c(round(nfeats*0.4), round(nfeats*0.4*0.6)), linear.output = FALSE)
    pred <- compute(m, subset(test, select = -state))
    p <- as.vector(pred$net.result)
    p_class <- ifelse(p>0.5, states[2], states[1])
    p_class <- factor(p_class, levels = states, labels = states)
    
  } else {
    
    stop('Invalid model type!')
    
  }
  
  return(list(predictions = p_class, observations = test$state, model = m))
  
}



measure.performance <- function(pred_class, obs_class) {
  
  print(paste0('Define ', toupper(states[2]), ' as Positive.'))
  
  perf <- confusionMatrix(pred_class, obs_class)

  if (perf[['positive']] == states[1]){
    tnr <- perf[[c('byClass','Sensitivity')]]
    tpr <- perf[[c('byClass','Specificity')]]
  } else if (perf[['positive']] == states[2]){
    tpr <- perf[[c('byClass','Sensitivity')]]
    tnr <- perf[[c('byClass','Specificity')]]
  }
  
  return(list(accuracy = perf[[c('overall','Accuracy')]], kappa = perf[[c('overall','Kappa')]], 
              sensitivity = tpr, specificity = tnr))
  
}



split.sample <- function(data, nfold = 10, mincount = 10){ 
 
  #set.seed(54)
  
  rawIdx <- list()
  count <- c()
  ndataPerFold <- c()
  
  for (si in 1:length(states)){ 
    
    rawIdx[[si]] <- c(1:nrow(data))[data$state == states[si]]
    count[si] <- length(rawIdx[[si]])
    
    ndataPerFold[si] <- ceiling(count[si]/nfold)
    
    if (ndataPerFold[si] <= ndataPerFold[si]*nfold - count[si]) {
      ndataPerFold[si] <- floor(count[si]/nfold)
    }
    
    if (count[si] > 1) {rawIdx[[si]] <- sample(rawIdx[[si]], count[si])} # shuffle idx
    
  }
  
  # enough samples?
  if (min(count) < mincount) {
    print(paste0('Class size less than ', mincount, '. Fail to split data.'))
    return(list())
  }
  
  trainIdx <- list()
  testIdx  <- list()
  
  for (foldi in 1:nfold){
    
    startPosition <- ndataPerFold * (foldi - 1) + 1
    
    if (foldi < nfold) {
      endPosition <- ndataPerFold * foldi
    } else {
      endPosition <- count
    }
    
    for (si in 1:length(states)){
      testPositions <- startPosition[si]:endPosition[si]
      if (si == 1){
        testIdx[[foldi]]  <- rawIdx[[si]][testPositions]
        trainIdx[[foldi]] <- rawIdx[[si]][-testPositions]
      } else {
        testIdx[[foldi]]  <- c(testIdx[[foldi]], rawIdx[[si]][testPositions])
        trainIdx[[foldi]] <- c(trainIdx[[foldi]], rawIdx[[si]][-testPositions])
      }
    }
    
  }
  
  return(list(trainIdx = trainIdx, testIdx = testIdx))
  
}



balance.class <- function(data, method = 'copy'){
  
  count <- c()
  for (si in 1:length(states)){
    eval(parse(text = paste0('data', si, ' <- subset(data, state == states[si])')))
    eval(parse(text = paste0('count[si] <- nrow(data', si, ')')))
  }
  
  # if one class is empty
  if (min(count) == 0){
    print('One of the classes is empty! Return NaN')
    return(NaN)
  }
  if (max(count) == min(count)) {return(data)}
  
  nCopy   <- floor(max(count)/min(count)) - 1
  nSelect <- max(count) %% min(count)
  
  if (method == 'copy'){
   
    data2copy <- eval(parse(text = paste0('data', which(count == min(count)))))
    
    if (nCopy > 0){
      for (copyi in 1:nCopy){
        if (copyi == 1){
          copy <- data2copy
        } else {
          copy <- rbind(copy, data2copy)
        }
      }
    }
    
    if (nSelect > 0){
      if (nCopy > 0){
        copy <- rbind(copy, data2copy[sample(1:min(count), nSelect),])
      } else {
        copy <- data2copy[sample(1:min(count), nSelect),]
      }
    } 
    
    newData <- rbind(data, copy)
    
  } else if (toupper(method) == 'SMOTE'){
    
    newData <- SMOTE(state~., data, perc.over = nCopy*100 + 100*(min(1,nSelect)), 
                     k = 5, perc.under = (1 + 1/(nCopy + min(1,nSelect))) * 100)
    
  } 
 
  return(newData)
  
}



normalize <- function(data, algorithm, pars = list()){
  
  # normalize data within each column
  # algorithm options are: 
  # - range
  # - z (pars: colMeans, colSds)
  
  # normalize each column (feature)
  labs <- data$state
  
  dataMat <- as.matrix(subset(data, select = -state))
  nObs    <- nrow(dataMat)
  
  if (algorithm == 'range'){
    
    minMat   <- matrix(rep(colMins(dataMat), nObs), nrow = nObs, byrow = TRUE)
    maxMat   <- matrix(rep(colMaxs(dataMat), nObs), nrow = nObs, byrow = TRUE)
    dataNorm <- (dataMat - minMat) / (maxMat - minMat)
    dataPars <- list(mins = colMins(dataMat), maxs = colMaxs(dataMat))
    
    
  } else if (algorithm == 'z'){
    
    if (length(pars) == 0) {
      meanMat  <- matrix(rep(colMeans(dataMat), nObs), nrow = nObs, byrow = TRUE)
      sdMat    <- matrix(rep(colSds(dataMat), nObs), nrow = nObs, byrow = TRUE)
    } else {
      meanMat  <- matrix(rep(pars$means, nObs), nrow = nObs, byrow = TRUE)
      sdMat    <- matrix(rep(pars$sds, nObs), nrow = nObs, byrow = TRUE)
    }
    dataNorm <- (dataMat - meanMat) / sdMat
    dataPars <- list(means = colMeans(dataMat), sds = colSds(dataMat))
    
  } else {
    
    return('Invalid algorithm!')
    
  }
  
  dataNorm <- as.data.frame(dataNorm)
  dataNorm$state <- labs
  
  return(list(dataNorm = dataNorm, dataPars = dataPars))
  
}



get.all.data <- function(subs, task, measures, feats, folders, probeLabOn, timeLabOn, nBack, states, states.train, normOn = TRUE){
  
  for (subi in 1:length(subs)){
    
    sub  <- subs[subi]
    temp <- get.data(sub, task, measures, feats, folders, probeLabOn, timeLabOn, nBack, states, states.train)
    print(paste0(nrow(temp), ' data have been loaded'))
    
    if (subi == 1) {
      rowCount  <- nrow(temp)
    } else {
      rowCount  <- rowCount + nrow(temp)
    }
    
    if (normOn) {
      temp <- normalize(temp, 'z')$dataNorm
    }
    
    if (subi == 1){
      data <- temp
    } else {
      data <- rbind(data, temp)
    }
    
  }  # loop over subs
  
  print(paste0('Overall row count is ', rowCount))
  return(data)
  
}



get.trialIdx <- function(sub, task, measures, feats, folders){
  
  for (feati in 1:length(feats)){
    
    feat   <- feats[feati]
    folder <- folders[feati]
    
    all <- readMat(paste0(f_measure_matfile, f_sep, folder, f_sep, sub, '.mat'))
    
    for (si in 1:length(states)){
      
      state <- states[si]
      eval(parse(text = paste0('temp <- all$', feat, '.', task, '.', state)))
      
      if (feati == 1){
        if (feat %in% c('alpha', 'theta')){
          temp <- cbind(temp[, 1:4], si - 1)  # trial idx + feats
        } else {
          temp <- cbind(temp[, 1:5], si - 1)
        }
      } else {
        if (feat %in% c('alpha', 'theta')){
          temp <- cbind(temp[, 3:4], si - 1)  # feats
        } else {
          temp <- cbind(temp[, 3:5], si - 1)
        }
      }
      
      if (si == 1){
        data <- temp 
      } else {
        data <- rbind(data, temp)
      }
      
    }
    
    if (feati == 1) {
      if (feat %in% c('alpha','theta')){
        colnames(data) <- c('ID.s1', 'ID.s2', paste0(measures[feati], '.', c('Base','StimOn')), 'state')
      } else {
        colnames(data) <- c('ID.s1', 'ID.s2', paste0(measures[feati], '.', c('Size','Time','Scale')), 'state')
        data[data[, 3] == 0, 1] <- NaN  # mark the detection failure of single trial ERP
      }
    } else {
      if (feat %in% c('alpha','theta')){
        colnames(data) <- c(paste0(measures[feati], '.', c('Base','StimOn')), 'state')
      } else {
        colnames(data) <- c(paste0(measures[feati], '.', c('Size','Time','Scale')), 'state')
        data[data[, 1] == 0, 1] <- NaN  # mark the detection failure of single trial ERP
      }
    }
    
    if (feati == 1){
      df <- data
    } else {
      df <- cbind(df, data[, -ncol(data)])
    }
    
  }
  
  df       <- as.data.frame(df)
  df$state <- factor(df$state, levels = c(1:length(states)) - 1, labels = states)
  df <- na.omit(df)  # remove feat detection failure trials
  df.idx <- df[,1:2] 
  
  return(df.idx)
}



permutate.data <- function(data){
  
  print(paste0('Shuffle labels of classes: ', paste0(toupper(states), collapse = ', ')))

  data <- data[sample.int(nrow(data)),]
  sizes <- summary(data$state)
  states <- names(sizes)  
  for (i in 1:(length(sizes))){
    state <- states[i]
    size  <- sizes[i]
    if (i == 1){
      data[1:size, 'state'] <- state
    } else {
      data[(1:size) + sum(sizes[1:(i-1)]), 'state'] <- state
    }
  }
  # note the state column type will be the same
  return(data)
}



get.data <- function(sub, task, measures, feats, folders, probeLabOn = FALSE, timeLabOn = FALSE, 
                     nBack = NaN, states, states.train = c(), printOn = FALSE){
  # states.train is the class for reading
  # states is the final class (relabelling: states.train -> states)
  
  if (length(states.train) == 0) {
    states.train <- states
  }
  
  print(paste0('Loading data ', sub, ' of task ', task, '...'))
  
  if (probeLabOn) {
    print('Label data based on self-reports.')
  } else if (timeLabOn) {
    print('Label data based on time course.')
  } else {
    print('Keep the original labels in the mat file.')
  }
  
  beh <- get.beh.data(sub, task, nBack)
  
  for (feati in 1:length(feats)){

    feat   <- feats[feati]
    folder <- folders[feati]
    ncolMea <- measureNcol[feati]
    if (printOn) {
      print(paste0('Loading marker ', feati, ' out of ', length(feats), ': ', measures[feati]))
    }

    all <- readMat(file.path(f_measure_matfile, folder, paste0(sub, '.mat')))
    
    for (si in 1:length(states)){
      
      state <- states[si]
     # if (!probeLabOn && !timeLabOn) {
     #    if (length(states.train) > 0) {
    #      eval(parse(text = paste0('temp <- all$', feat, '.', states.train[si])))
     #   } else {
    #      eval(parse(text = paste0('temp <- all$', feat, '.', state)))
     #   }
    #  } else {
      if (probeLabOn) {
        urevents <- beh[beh$state == state, 'urevent']
      } else if (timeLabOn) {
        urevents <- beh[beh$fatigue == states.train[si], 'urevent']
      } else {
        urevents <- beh[beh$block == states.train[si], 'urevent']
      }
      
      if (length(urevents) > 0) {
        for (li in 1:length(all)) {  # loop over all data
          tp <- all[[li]]
          tp <- tp[tp[,1] %in% urevents, ]  # filter based on urevents
          if (li == 1) {
            temp <- tp
          } else {
            temp <- rbind(temp, tp)
          }
        }
      } 
     # } 
      
      # select feats, add lab. Note cbind() wrongly handle one-row data
      if ((!probeLabOn && !timeLabOn) || (length(urevents)>0 && nrow(temp) > 0)) { 
        temp2 <- temp[,1] # save the index
        if (feat %in% c('alpha', 'theta')){
          if (nrow(temp) > 1) {temp <- cbind(temp[, 2:3], si - 1)} else {temp <- matrix(c(temp[2:3], si-1), nrow = 1, byrow = TRUE)}
        } else if (feat %in% c('p1', 'n1', 'p3')){
          if (nrow(temp) > 1) {temp <- cbind(temp[, 2:4], si - 1)} else {temp <- matrix(c(temp[2:4], si-1), nrow = 1, byrow = TRUE)}
        } else {
          if (nrow(temp) > 1) {temp <- cbind(temp[, 2:ncolMea], si - 1)} else {temp <- matrix(c(temp[2:ncolMea], si-1), nrow = 1, byrow = TRUE)}
        } 
      } else { # no such condition?
        temp2 <- c()
        if (feat %in% c('alpha', 'theta')){
          temp <- matrix(c(NaN, NaN, si - 1), byrow = TRUE, nrow = 1)
        } else if (feat %in% c('p1', 'n1', 'p3')){
          temp <- matrix(c(NaN, NaN, NaN, si - 1), byrow = TRUE, nrow = 1)
        } else {
          temp <- matrix(c(rep(NaN, ncolMea-1), si - 1), byrow = TRUE, nrow = 1)
        } 
        temp <- na.omit(temp)
      }
      
      if (si == 1){
        data  <- temp 
        urevs <- temp2
      } else {
        if (nrow(data) > 0){
          if (nrow(temp) > 0){
            data <- rbind(data, temp)
            urevs <- c(urevs, temp2)
          }  # else, no act
        } else {
          data <- temp
          urevs <- temp2
        }
      }
      
    }  # loop over states
    
    # feature combination
    if (nrow(data) > 0){
      if (feat %in% c('alpha','theta')) {
        colnames(data) <- c(paste0(measures[feati], '.', c('pre','post')), 'state')
      } else if (grepl('ica', feat) || grepl('pca', feat)) {
        if (ncolMea == 3) {
          colnames(data) <- c(paste0(measures[feati], '.', c('pre','post')), 'state')
        } else {
          colnames(data) <- c(paste0(measures[feati], '.pnt', 1:(ncol(data)-1)), 'state')
        }
      } else if (feat %in% c('p1', 'n1', 'p3')){
        colnames(data) <- c(paste0(measures[feati], '.', c('Size','Time','Scale')), 'state')
        data[data[, 1] == 0, 1] <- NaN  # mark the detection failure of single trial ERP
      }
      
      # rearrange data by the urev order
      urorder <- order(urevs)
      urevs <- urevs[urorder]
      data  <- data[urorder,]
      
      if (feati == 1){
        df <- data
        urevs1 <- urevs
      } else {
        if (all(urevs1 == urevs)) {
          df <- cbind(df, data[, -ncol(data)])
        } else {
          stop('Feature urevent inconsisdent!')
        }
      }
    } else { # if data is empty for this feature, same for the other feats
      df <- data
    }

  }  # loop over feats
  
  df <- as.data.frame(df)
  
  # process labels
  df$state <- factor(df$state, levels = c(1:length(states)) - 1, labels = states)
  
  # remove na
  df <- na.omit(df)  # remove feat detection failure trials
  
  print('Successful.')
  return(df)
  
}



simulate.data <- function(sub, task, measures, measureTypes, hIncreaseOn){
  
  load('trialCount.rdata')
  centers <- c(0, 1)
  sharp   <- 1
  
  eval(parse(text = paste0('trialCount <- trialCount.', task)))
  
  for (si in 1:length(states)){
    state <- states[si]
    
    # get class size
    size <- trialCount[sub, state]
    
    # simlute data
    for (mi in 1:length(measures)){
      
      measure     <- measures[mi]
      measureType <- measureTypes[mi]
      increaseOn  <- hIncreaseOn[mi]
      
      if (state == 'ot'){
        if (increaseOn) {
          center <- min(centers)
        } else {
          center <- max(centers)
        }
      } else {
        if (increaseOn){
          center <- max(centers)
        } else {
          center <- min(centers)
        }
      }
      
      if (toupper(measureType) == 'ERP'){
        temp <- cbind(rnorm(size, mean = center, sd = sharp),
                      rnorm(size, mean = mean(centers), sd = sharp),
                      rnorm(size, mean = mean(centers), sd = sharp))
        colnames(temp) <- paste0(measure, '.', c('Size','Time','Scale'))
      } else {
        temp <- cbind(rnorm(size, mean = center, sd = sharp),
                      rnorm(size, mean = center, sd = sharp))
        colnames(temp) <- paste0(measure, '.', c('Base','StimOn'))
      }
      
      if (mi == 1){
        data <- temp
      } else {
        data <- cbind(data, temp)
      }
 
    }
    
    data <- cbind(data, state = si-1)
    if (si ==  1){
      df <- data
    } else {
      df <- rbind(df, data)
    }
    
  }

  df       <- as.data.frame(df)
  df$state <- factor(df$state, levels = c(1:length(states)) - 1, labels = states)

  return(df)
  
}



simulate.data.random <- function(sub, task, measures, measureTypes, equalSizeOn = FALSE){
  
  center <- 0
  sharp  <- 1
  size.default <- 100
  
  if(is.logical(equalSizeOn) && !equalSizeOn) {
    load('trialCount.rdata')
    eval(parse(text = paste0('trialCount <- trialCount.', task)))
  }
  
  for (si in 1:length(states)){
    state <- states[si]
    
    # get class size
    if (is.logical(equalSizeOn) && !equalSizeOn){
      size <- trialCount[sub, state]
    } else if (is.logical(equalSizeOn) && equalSizeOn) {
      size <- size.default
    } else if (is.numeric(equalSizeOn)){
      size <- equalSizeOn
    } else {
      size <- NaN
    }
    
    
    # simlute data
    for (mi in 1:length(measures)){
      
      measure     <- measures[mi]
      measureType <- measureTypes[mi]
   
      if (toupper(measureType) == 'ERP'){
        temp <- cbind(rnorm(size, mean = center, sd = sharp),
                      rnorm(size, mean = center, sd = sharp),
                      rnorm(size, mean = center, sd = sharp))
        colnames(temp) <- paste0(measure, '.', c('Size','Time','Scale'))
      } else {
        temp <- cbind(rnorm(size, mean = center, sd = sharp),
                      rnorm(size, mean = center, sd = sharp))
        colnames(temp) <- paste0(measure, '.', c('Base','StimOn'))
      }
      
      if (mi == 1){
        data <- temp
      } else {
        data <- cbind(data, temp)
      }
      
    }
    
    data <- cbind(data, state = si-1)
    if (si ==  1){
      df <- data
    } else {
      df <- rbind(df, data)
    }
    
  }
  
  df       <- as.data.frame(df)
  df$state <- factor(df$state, levels = c(1:length(states)) - 1, labels = states)
  
  return(df)
  
}



filter.feats <- function(keyCombos = list(), feats2load = c()) {
  
  source(file.path(f_code, 'pars.r'))
  
  if (length(keyCombos) > 0) {
    for (fi in 1:length(keyCombos)) {
      
      featType <- keyCombos[[fi]][1]
      measureType <- keyCombos[[fi]][2]
      
      if (featType != '') {
        tempIdx <- feats == featType
      } else {
        tempIdx <- rep(TRUE, length(feats))
      }
      
      if (measureType != '') {
        tempIdx <- tempIdx & measureTypes == measureType
      } 
      
      # combine idx
      if (fi == 1){
        idx <- tempIdx
      } else {
        idx <- idx | tempIdx
      }
      
    }  # loop over keyCombos
  } else {
    idx <- rep(TRUE, length(feats))
  }
  
  if (length(feats2load) > 0) {
    idx <- idx & measures %in% feats2load
  }
  
  feats <- feats[idx]
  measures <- measures[idx]
  measureNames <- measureNames[idx]
  folders <- folders[idx]
  measureTypes <- measureTypes[idx]
  measureColors <- measureColors[idx]
  measureShapes <- measureShapes[idx]
  hIncreaseOn <- hIncreaseOn[idx]
  measureNcol <- measureNcol[idx]
  
  save(file = 'temp_pars.rdata',feats, measures, measureNames, folders, 
       measureTypes, measureColors, measureShapes, hIncreaseOn, measureNcol)
  load('temp_pars.rdata', .GlobalEnv)
}
featTesting <- function(subs, tasks, mtype, validType, parMat,
                        taskMerge, useParMat, permutateOn = FALSE,
                        gridSearchOn = FALSE, searchBy = '', 
                        parallelProcId = 0, normalizeOn, groupOn = FALSE,
                        lopocvOn = FALSE) {
# can only be performed on successfully modelled data
  
  if (mtype == 'svm' && gridSearchOn){
    cPowerVec <- -5:15
    gPowerVec <- -15:3
    
    print('Initialize par matrix.')
    if (groupOn && !lopocvOn) {
      print('Group modelling is ON. Training data of ALL the participants.')
      parMat <- rep(0, length(feats) * 2)
      names(parMat) <- paste0(rep(measures, each = 2) , '_',  c('c', 'gamma'))
    } else {
      if (groupOn && lopocvOn) {
        print('Leave-one-participant-out cross-validation is On. Training data of (N-1) participants.')
      } else { 
          print('Group modelling is OFF. Training data in EACH participant.')
      }
      parMat <- matrix(0, length(subs), length(feats) * 2)
      colnames(parMat) <- paste0(rep(measures, each = 2) , '_',  c('c', 'gamma'))
    }

  }
  
  if (length(tasks) == 1 || taskMerge) {
    
    print('Fit model on task-combined data.')
    
    print('Initialze accuracy matrix.')
    if (groupOn && !lopocvOn) {
      accPerFeat <- rep(0, length(feats))
      names(accPerFeat)  <- measures
    } else {
      accPerFeat <- matrix(0, length(subs), length(feats))
      colnames(accPerFeat) <- measures
      rownames(accPerFeat) <- subs
    }
    
    if (!groupOn) {
      
      # monitor progress
      n  <- length(subs) * length(feats)
      ti <- 0
      if (!is.numeric(parallelProcId)) {
        pb0 <- tkProgressBar(title = 'Testing each feature', min = 0, max = n, width = 300)
      } else {
        pb0 <- tkProgressBar(title = paste0('Testing each feature', parallelProcId), min = 0, max = n, width = 300)
      }
      
      for (subi in 1:length(subs)){
        
        sub <- subs[subi]
        if (mtype == 'svm' && !gridSearchOn && useParMat){
          if (ncol(parMat) > 2){
            parLists <- list()
            for (taski in 1:length(tasks)) {
              task <- tasks[taski]
              parLists[[taski]] <- list(c = parMat[subi, paste0(task, '_c')], gamma = parMat[subi, paste0(task, '_c')])
            }
            
          } else {
            parList <- list(c = parMat[subi, 'c'], gamma = parMat[subi, 'gamma'])
          }
        } else {
          print('Grid search is off and no specified pars. Use default pars to model!')
          parList <- list()
        }
        tic(paste('Testing each feature on sub', sub))
        
        for (feati in 1:length(feats)){
          
          print(paste0('Test marker ', measureNames[feati]))
          for (taski in 1:length(tasks)){
            task <- tasks[taski]
            
            # load
            temp <- get.data(sub, task, measures[feati], feats[feati], folders[feati], 
                             probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack, 
                             states = states, states.train = states.train)
            
            # normalize: within tasks
            if (normalizeOn) {
              print('Normlize loaded feature within tasks.')
              temp <- normalize(temp, 'z')$dataNorm
              if (taski == 1){
                data <- temp  # intialize
              } else {
                print('Combine across tasks.')
                data <- rbind(data, temp)
              }
            }
            
          }  # over tasks
          
          # shuffle labels
          if (permutateOn){
            data <- permutate.data(data)
          }
          
          # get parList, fit models
          if (mtype == 'svm' && !gridSearchOn && useParMat){ # grid search off, parMat specified
            if (ncol(parMat) > 2){
              perfs <- list()
              for (taski in 1:length(tasks)) {
                if (validType == 'loo') {
                  perfs[[taski]] <- leave.one.out(data, mtype, parList = parLists[[taski]], balanced = 'copy')
                } else {
                  perfs[[taski]] <- cross.validation(data, mtype, parList = parLists[[taski]], balanced = 'copy')
                  
                }
              }
              
            } else {
              if (validType == 'loo'){
                perf <- leave.one.out(data, mtype, parList = parList, balanced = 'copy')
              } else {
                perf <- cross.validation(data, mtype, parList = parList, balanced = 'copy')
              }
            }
            
          } else if (mtype == 'svm' && gridSearchOn){ # grid search on
            
            mat <- grid.search.svm(data, 2^cPowerVec, 2^gPowerVec, validType = 'cv', searchBy = searchBy)
            # fail to do grid search?
            if (!is.matrix(mat) && is.nan(mat)){
              print('Fail to do grid search. Use default parameters to fit model!')
              parMat[subi, c(1:2) + (feati - 1)*2] <- NaN
              parList <- list()
            } else {
              pos <- which.localMax(mat, gPowerVec, cPowerVec)
              parMat[subi, c(1:2) + (feati - 1)*2] <- 2^pos
              parList <- list(c = parMat[subi, paste0(measures[feati], '_c')], gamma = parMat[subi, paste0(measures[feati], '_gamma')])
            }
            
            if (validType == 'loo'){
              perf <- leave.one.out(data, mtype, parList = parList) # [default] balanced = 'copy' for the training split
            } else {
              perf <- cross.validation(data, mtype, parList = parList) # [default] balanced = 'copy' for the training split
            }
            
          } else { # grid search off, no parMat specified
            if (validType == 'loo'){
              perf <- leave.one.out(data, mtype, balanced = 'copy')
            } else {
              perf <- cross.validation(data, mtype, balanced = 'copy')
            }
          }  # end of model fitting
          
          if (mtype == 'svm' && !gridSearchOn && useParMat && ncol(parMat) > 2){
            acc_vec <- c()
            for (taski in 1:length(tasks)) {
              acc_vec[taski] <- max(perfs[[taski]]$accuracy)
            }
            accPerFeat[subi, feati] <- max(acc_vec)
          } else {
            accPerFeat[subi, feati] <- perf$accuracy
          }
          
          if (!is.numeric(parallelProcId)) {
            f_temp <- 'temp_featTesting.rdata'
          } else {
            f_temp <- paste0('temp_featTesting', parallelProcId, '.rdata')
          }
          
          print(paste0('Autosave in ', f_temp))
          if (mtype == 'svm' && gridSearchOn){
            save(file = f_temp, accPerFeat, parMat)
          } else {
            save(file = f_temp, accPerFeat)
          }
          
          ti <- ti + 1
          setTkProgressBar(pb0, ti, label = paste(round(ti/n * 100, 1), "% done"))
          
        }  # loop over feats
        
        toc()
        
      }  # loop over subs
      close(pb0)
    
    } else if (groupOn && !lopocvOn){
      
      for (feati in 1:length(feats)){
        print(paste0('Test marker ', measureNames[feati]))
        tic(paste0('Test marker ', measureNames[feati], ' on ALL participants'))
        
        for (subi in 1:length(subs)){
          
          sub <- subs[subi]
          
          for (taski in 1:length(tasks)){
            task <- tasks[taski]
            
            # load
            temp <- get.data(sub, task, measures[feati], feats[feati], folders[feati], 
                             probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack, 
                             states = states, states.train = states.train)
            
            # normalize: within tasks
            if (normalizeOn) {
              print('Normlize loaded feature per task/individual.')
              temp <- normalize(temp, 'z')$dataNorm
              if (taski == 1){
                data.tp <- temp  # intialize
              } else {
                print('Combine across tasks.')
                data.tp <- rbind(data.tp, temp)
              }
            }
            
          }  # over tasks
          print(paste0(nrow(data.tp), ' data has been loaded.'))
          
          if (subi == 1) {
            
            data <- data.tp
          } else {
            print('Combine across individuals.')
            data <- rbind(data, data.tp)
            print(paste0('Overall data count is ', nrow(data)))
          }
          
        }  # loop over subs
          
        # shuffle labels
        if (permutateOn){
          data <- permutate.data(data)
        }
          
        # get parList, fit models
        # progress monitor
        print(paste0('Progress ', feati, ' out of ', length(feats)))
        
        if (mtype == 'svm' && gridSearchOn){ # grid search on
          
          print(paste0('Grid search for best SVM parameters base on ', toupper(searchBy)))
          mat <- grid.search.svm(data, 2^cPowerVec, 2^gPowerVec, validType = 'cv', searchBy = searchBy)
          # fail to do grid search?
          if (!is.matrix(mat) && is.nan(mat)){
            print('Fail to do grid search. Use default parameters to fit model!')
            parMat[c(1:2) + (feati - 1)*2] <- NaN
            parList <- list()
          } else {
            pos <- which.localMax(mat, gPowerVec, cPowerVec)
            parMat[c(1:2) + (feati - 1)*2] <- 2^pos
            parList <- list(c = parMat[paste0(measures[feati], '_c')], gamma = parMat[paste0(measures[feati], '_gamma')])
          }
          
          if (validType == 'loo'){
            perf <- leave.one.out(data, mtype, parList = parList) # [default] balanced = 'copy' for the training split
          } else {
            perf <- cross.validation(data, mtype, parList = parList) # [default] balanced = 'copy' for the training split
          }
          
        } else { # grid search off, no parMat specified
          print('Grid search is OFF and no pars specified.')
          print('Use DEFAULT pars to fit model.')
          
          if (validType == 'loo'){
            perf <- leave.one.out(data, mtype, balanced = 'copy')
          } else {
            perf <- cross.validation(data, mtype, balanced = 'copy')
          }
        }  # end of model fitting
        
        # register
        accPerFeat[feati] <- perf$accuracy
        
        # auto temp save
        if (!is.numeric(parallelProcId)) {
          f_temp <- 'temp_featTesting.rdata'
        } else {
          f_temp <- paste0('temp_featTesting', parallelProcId, '.rdata')
        }
        
        print(paste0('Autosave in ', f_temp))
        if (mtype == 'svm' && gridSearchOn){
          save(file = f_temp, accPerFeat, parMat)
        } else {
          save(file = f_temp, accPerFeat)
        }
        
        toc()
        
      } # loop over feats
      
    # end of groupOn   
    } else if (groupOn && lopocvOn){
      
      for (feati in 1:length(feats)){
        print(paste0('Test marker ', measureNames[feati]))
        
        for (subi in 1:length(subs)){
          
          sub <- subs[subi]
          tic(paste0('Test marker ', measureNames[feati], ' on participant ', sub))
          
          for (taski in 1:length(tasks)){
            task <- tasks[taski]
            
            # load
            temp <- get.all.data(setdiff(subs,sub), task, measures[feati], feats[feati], folders[feati], 
                                 probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack, 
                                 states = states, states.train = states.train, normOn = normalizeOn)
            temp2 <- get.all.data(sub, task, measures[feati], feats[feati], folders[feati], 
                                  probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack, 
                                  states = states, states.train = states.train, normOn = normalizeOn)
            if (normalizeOn) {print('Normlize loaded feature per task/individual.')}
            
            if (taski == 1){
              data.train <- temp
              data.test  <- temp2
            } else {
              print('Combine across tasks.')
              data.train <- rbind(data.train, temp)
              data.test  <- rbind(data.test, temp2)
            }
            
          }  # over tasks
          
          # shuffle labels
          if (permutateOn){
            data.train <- permutate.data(data.train)
            # ??? shuffle data.test too?
          }
          
          # balance training dataset
          data.train.balan <- balance.class(data.train, 'copy')
          
          # get parList, fit models
          # progress monitor
          print(paste0('Progress ', round(((feati-1)*length(feats)+subi)/length(subs)/length(feats)*100, 2), '%'))
          
          if (mtype == 'svm' && gridSearchOn){ # grid search on
            
            print(paste0('Grid search for best SVM parameters base on ', toupper(searchBy)))
            mat <- grid.search.svm(data, 2^cPowerVec, 2^gPowerVec, validType = 'cv', searchBy = searchBy)
            # fail to do grid search?
            if (!is.matrix(mat) && is.nan(mat)){
              print('Fail to do grid search. Use default parameters to fit model!')
              parMat[subi, c(1:2) + (feati - 1)*2] <- NaN
              parList <- list()
            } else {
              pos <- which.localMax(mat, gPowerVec, cPowerVec)
              parMat[subi, c(1:2) + (feati - 1)*2] <- 2^pos
              parList <- list(c     = parMat[subi, paste0(measures[feati], '_c')], 
                              gamma = parMat[subi, paste0(measures[feati], '_gamma')])
            }
            
            
            m <- feat.modeling(data.train.balan, data.test, mtype, parList = parList)
            perf <- measure.performance(m$predictions, m$observations)
            
          } else { # grid search off, no parMat specified
            print('Grid search is OFF and no pars specified.')
            print('Use DEFAULT pars to fit model.')
            
            m <- feat.modeling(data.train.balan, data.test, mtype)
          }  # end of model fitting
          
          # register
          perf <- measure.performance(m$predictions, m$observations)
          accPerFeat[subi, feati] <- perf$accuracy
          
          # auto temp save
          if (!is.numeric(parallelProcId)) {
            f_temp <- 'temp_featTesting.rdata'
          } else {
            f_temp <- paste0('temp_featTesting', parallelProcId, '.rdata')
          }
          
          print(paste0('Autosave in ', f_temp))
          if (mtype == 'svm' && gridSearchOn){
            save(file = f_temp, accPerFeat, parMat)
          } else {
            save(file = f_temp, accPerFeat)
          }
          
          toc()
          
        }  # loop over subs
        
      } # loop over feats
      
    } # end of groupOn && lopocvOn
    
    
  } else { # task-separate
    
    print('Fit model on each task.')
    
    accPerFeat <- matrix(0, length(subs), length(feats))
    colnames(accPerFeat) <- measures
    rownames(accPerFeat) <- subs

    for (taski in 1:length(tasks)){
      eval(parse(text = paste0('accPerFeat.', tasks[taski], ' <- accPerFeat')))
    }
    
    for (taski in 1:length(tasks)){
      
      task <- tasks[taski]
      
      n <- length(subs) * length(feats)
      ti <- 0
      pb0 <- tkProgressBar(title = paste('Testing each feature of', task), min = 0, max = n, width = 300)
      
      for (subi in 1:length(subs)){
        
        sub <- subs[subi]
        tic(paste('Testing each feature of', task, 'on subject', sub))
        
        for (feati in 1:length(feats)){
          
          # load
          data <- get.data(sub, task, measures[feati], feats[feati], folders[feati], 
                           probeLabOn = probeLabOn.train, timeLabOn = timeLabOn.train, nBack, 
                           states = states, states.train = states.train)
          
          # normalize
          if (normalizeOn) {
            data <- normalize(data, 'z')$dataNorm
          }
          
          # permutate
          if (permutateOn) {data <- permutate.data(data)}
          
          # specified pars?
          if(useParMat) {
            parList <- list(c = parMat[subi, paste0(task, '_c')], gamma = parMat[subi, paste0(task, '_gamma')])
          } else {
            parList = list()
          }
          
          if (validType == 'loo'){
            perf <- leave.one.out(data, mtype, parList = parList, balanced = 'copy' )
          } else {
            perf <- cross.validation(data, mtype, parList = parList, balanced = 'copy')
          }
          
          eval(parse(text = paste0('accPerFeat.', task, '[subi, feati] <- perf$accuracy')))
          
          if (parallelProcId == 0) {
            f_temp = 'temp_featTesting.rdata'
          } else {
            f_temp = paste0('temp_featTesting', parallelProcId, '.rdata')
          }
          save(file = f_temp, list = paste0('accPerFeat.', tasks))
          print(paste0('Autosave in ', f_temp))
          
          ti <- ti + 1
          setTkProgressBar(pb0, ti, label = paste(round(ti/n * 100, 1), "% done"))
        }
        toc()
      }
      
      close(pb0)
    }
  }  # task-separate over
  
  return(f_temp)
  
}
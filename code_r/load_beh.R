library(dplyr)
library(R.matlab)
library(stringr)

# specify the folder of the csvs
f_raw <- 'C:\\TOPIC_mind wandering\\3data2\\data_share\\raw'  

# e.g., to get data of three preceding trials in SART of dataset 301
# df <- get.beh.data(301, 'sart', 3)


get.beh.data <- function(sub, task, nBack = NaN) {
  
 
  if (task == 'sart') {
    data <- get.data.sart(sub)
  } else if (task == 'vs') {
    data <- get.data.vs(sub)
  } else {
    stop('Invalid task.')
  }
  
  if (!is.nan(nBack)) {
    print(paste0('Get data from preceding ', nBack, ' trials.'))
    data <- get.data.nbk(data, nBack)
  } else {
    print('Get full behavioral data.')
    }
  data <- add.labs(data)
  
}



add.labs <- function(data, nsplit = 2) {

  print(paste0('Split data into ', nsplit, ' by time course'))
  nblock <- max(data$blocks.thisN, na.rm = TRUE) + 1
  split  <- round(nblock/nsplit)
  data$fatigue <- NaN
  data$fatigue <- ifelse(data$blocks.thisN < split, 'nf', 
                         ifelse(data$blocks.thisN >= nblock-split, 'f', 'nfa'))
  
  data$state <- NaN
  #data <- subset(data, pr.resp != 0)  # filter out uncertain responses 
  data$state <- ifelse(data$pr.resp > 0, 'ot', ifelse(data$pr.resp < 0, 'mw', 'uc'))  # 'uc' - uncertain
  return(data)
}



get.data.nbk <- function(data, n) {
  
  data <- subset(data, dis2pr < n)
  return(data)
  
}



get.data.sart <- function(sub, ureventOn = FALSE) {
  

  # pars
  cols2analyze <- c('number', 'type', 'trigger', 'blocks.thisN', 'trials.thisN', 
                    'resp.keys', 'resp.corr', 'resp.rt', 
                    'rating.response', 'rating.rt')
  
  data.raw <- read.csv(file.path(f_raw, paste0(sub, '_sart.csv')), header = TRUE, stringsAsFactors = FALSE)
  
  # register sub info
  #register.subinfo(sub, data.raw)
  
  
  # filter out block_prac
  rmRowId   <- which(is.na(data.raw$rating.response) & is.na(data.raw$resp.corr))
  data <- data.raw[-rmRowId,]
  
  # gather columns
  data <- data[,cols2analyze]
  
  # add dis2pr
  data <- add.dis2pr(data)
  
  # add urevent
  #if (ureventOn) {data <- add.urevent(data, sub)}
  
  return(data)
}



get.prac.beh <- function(sub, task) {
  
  data.raw <- read.csv(file.path(f_raw, paste0(sub, '_', task, '.csv')), header = TRUE, stringsAsFactors = FALSE)
  if (task == 'vs') {
    data <- subset(data.raw, !is.na(blocks_prac.thisRepN))
  } else if (task == 'sart') {
    data <- subset(data.raw, !is.na(prac_blocks.thisRepN))
  }
  
  return(data)
  
}



get.data.vs <- function(sub, ureventOn = FALSE) {
  

  # pars
  cols2analyze <- c('block', 'trigger', 'target', 'blocks.thisRepN', 'blocks.thisN', 'trials.thisN', 
                    'nTri', 'nSqr', 'nPen', 'nHex', 'resp.keys', 'resp.corr', 'resp.rt', 
                    'rating.response', 'rating.rt')
  
  data.raw <- read.csv(file.path(f_raw, paste0(sub, '_vs.csv')), header = TRUE, stringsAsFactors = FALSE)
  
  # register sub info
  #register.subinfo(sub, data.raw)
  
  
  # filter out block_prac
  startRowId <- min(which(data.raw$blocks.thisRepN == 0))
  data <- data.raw[startRowId:nrow(data.raw),]
  
  # gather columns (some modification after pilots)
  if ('group' %in% colnames(data)) {cols2analyze <- c(cols2analyze, 'group')}
  if (!'nHex' %in% colnames(data)) {cols2analyze <- cols2analyze[cols2analyze != 'nHex']}
  data <- data[,cols2analyze]
  
  # add dis2pr
  data <- add.dis2pr(data)
  
  # add urevents
  #if (ureventOn) {data <- add.urevent(data, sub)}
  
  return(data)
}






add.dis2pr <- function(data) {
  data$dis2pr <- NaN
  data$pr.resp <- NaN
  data$pr.rt <- NaN
  
  blocks <- as.numeric(na.omit(unique(data$blocks.thisN)))
  probes <- na.omit(data[,c('rating.response','rating.rt')])
  
  # check
  if (length(blocks) != nrow(probes)) {
    stop('Unequal count of block and probe!')
  }
  
  for (bi in 1:length(blocks)) {
    block <- blocks[bi]
    rows <- data$blocks.thisN == block
    rows[is.na(rows)] <- FALSE
    trials <- data[rows,]$trials.thisN
    data[rows,]$dis2pr <- length(trials) - trials - 1
    data[rows, c('pr.resp', 'pr.rt')]  <- probes[bi,]
  }
  
  return(data)
}









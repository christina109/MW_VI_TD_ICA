library(dplyr)
library(ggplot2)
library(ggpubr)
library(R.matlab)
library(stringr)

# path
f_pars   <- file.path(f_main, 'pars.rdata')
f_raw    <- file.path(f_main, 'raw')

# pars
load(f_pars)

analyze.beh <- function(subs, ureventOn = FALSE, individualOn = TRUE, savePlotOn = TRUE){
  
  # intialize
  plist <- list()
  glist <- list()
  rlist <- list()
  alist <- list()
  
  # output name
  if (length(subs) > 1) {
    subchar <- paste0('_', min(subs), '_', max(subs))
  } else {
    subchar <- paste0('_', subs)
  }
  
  for (si in 1:length(subs)){
    sub  <- subs[si]
    
    # process decision
    sartOn <- FALSE
    if (file_test('-f', file.path(f_raw, paste0(sub, '_sart.csv')))) {
      sartOn <- TRUE
    } 
    
    data <- get.data.vs(sub, ureventOn)
    if (sartOn) {
      data.sart <- get.data.sart(sub, ureventOn)
    }
    
    # get line id
    probeId.record <- which(!is.na(data$rating.response))
    probeId.true   <- probeId.record - 1
    if (sartOn) {
      probeId.record.sart <- which(!is.na(data.sart$rating.response))
    }
    
    # combine data
    data.rating <- data.frame(cond = data[probeId.true, 'block'], rating = data[probeId.record, 'rating.response'])
    data.rating$block <- 1:nrow(data.rating)
    if (sartOn) {
      data.rating.sart <- data.frame(cond = 'sart', rating = data.sart[probeId.record.sart, 'rating.response'])
      data.rating.sart$block <- 1:nrow(data.rating.sart)
      data.rating <- rbind(data.rating, data.rating.sart)
    }
    
    # condition stat
    data.rating <- group_by(data.rating, cond)
    temp <- summarise(data.rating, mean.score = round(mean(rating), 2))
    for (condi in 1:nrow(temp)) {
      m <- lm(rating~block, subset(data.rating, cond == temp$cond[condi]))
      temp[condi, 'r2'] <- round(summary(m)$r.squared, 3)
      temp[condi, 'p']  <- summary(m)$coefficients[2,4]
      temp[condi, 'sig'] <- ifelse(temp[condi, 'p'] <= 0.001, '***', ifelse(temp[condi, 'p'] <= 0.01, "**", ifelse(temp[condi, 'p'] <= 0.05, "*", '')))
    }
    
    # rt data
    beh <- data[, c('block', 'target', 'resp.corr', 'resp.rt')]
    beh <- beh[1:(nrow(beh)-1),]
    beh[beh$block == 'ioa', 'target'] <- 'none'
    if (sartOn) {
      beh.sart <- data.frame(block = 'sart', target = data.sart$type,  data.sart[, c('resp.corr', 'resp.rt')])
      beh.sart <- beh.sart[1:(nrow(beh.sart)-1),]
      beh <- rbind(beh, beh.sart)
    }
    beh$target <- factor(beh$target, levels = c('triangle', 'square', 'pentagon', 'hexagon', 'none', 'nt', 'target'), 
                         labels = c('Triangle', 'Square', 'Pentagon', 'Hexagon', 'None', 'GO', 'NOGO'))
    
    # add to group data
    if (si == 1) {
      data.rating.all <- data.frame(data.rating, sub = sub)
      beh.all <- data.frame(beh, sub = sub)
    } else {
      data.rating.all <- rbind(data.rating.all, data.frame(data.rating, sub = sub) )
      beh.all <- rbind(beh.all, data.frame(beh, sub = sub))
    }
    
    # plot: rating bar
    p <- ggplot(data.rating, aes(x = factor(rating, levels = -5:5), fill = cond)) 
    p <- p + geom_bar(position = 'dodge')
    p <- p + scale_x_discrete(drop = FALSE)
    p <- p + labs(title = paste0('Data ', sub), x = 'Rating', y = 'Count', caption = 'Number in brackets indicates mean score')
    p <- p + scale_fill_manual('Block', labels = paste0(toupper(levels(data.rating$cond)), '(',temp$mean.score, ')'), 
                               values = as.character(cond.colors[which(cond.colors[,'cond']== levels(data.rating$cond)), 'color']))
    p <- p + theme(axis.text.x = element_text(color = rating.colors))
    
    # save to a list
    plist[[si]] <- p
    
    # plot: rating trend
    p <- ggplot(data.rating, aes(x = block, y = rating, col = cond))
    p <- p + geom_point()
    p <- p + geom_smooth(method = 'lm', se = FALSE) 
    p <- p + scale_y_continuous(limits = c(-5,5), breaks = -5:5)
    p <- p + labs(title = paste0('Data ', sub), x = 'Block sequence', y = 'Rating', caption = as.expression(bquote('Number in brackets indicates' ~ italic(R)^2)))
    p <- p + scale_color_manual('Block', labels = paste0(toupper(levels(data.rating$cond)), '(',temp$r2, temp$sig, ')'), 
                               values = as.character(cond.colors[which(cond.colors[,'cond']== levels(data.rating$cond)), 'color']))
    p <- p + theme(axis.text.y = element_text(color = rating.colors))
    
    # save to a list
    glist[[si]] <- p
    
    # plot: rt distribution
    p <- ggplot(beh[beh$resp.corr == 1 & beh$target != 'NOGO',], aes(x = target, y = resp.rt, col = target))
    p <- p + geom_point(alpha = 0.5)
    p <- p + stat_summary(fun.data = "mean_cl_normal", geom = "crossbar", width = 0.2) 
    p <- p + labs(title = paste0('Data ', sub), x = '', y = 'Mean RT (s)', caption = 'EOA includes the last three conditions')
    rlist[[si]] <- p
    
    # plot: acc distribution
    p <- ggplot(beh, aes(x = target, y = resp.corr, col = target))
    p <- p + geom_point(alpha = 0.5)
    p <- p + stat_summary(fun.data = "mean_cl_normal", geom = "crossbar", width = 0.2) 
    p <- p + labs(title = paste0('Data ', sub), x = '', y = 'Accuracy', caption = 'EOA includes the last three conditions')
    alist[[si]] <- p
    
  }  # loop over subs
  
  # stat: rating
  print('Preparing data to plot...')
  data.rating.all <- group_by(data.rating.all, sub, cond)
  data.mean <- summarise(data.rating.all, mean.score = round(mean(rating),2))
  temp <- data.frame(cond = levels(data.rating.all$cond))
  for (condi in 1:nrow(temp)){
    m <- lm(rating~block, subset(data.rating.all, cond == temp$cond[condi]))
    temp[condi, 'r2'] <- format(summary(m)$r.squared, digits = 3)
    temp[condi, 'p']  <- summary(m)$coefficients[2,4]
    temp[condi, 'sig'] <- ifelse(temp[condi, 'p'] <= 0.001, '***', ifelse(temp[condi, 'p'] <= 0.01, "**", ifelse(temp[condi, 'p'] <= 0.05, "*", '')))
  }
  
  # stat: beh
  beh.all  <- group_by(beh.all, sub, target)
  beh.stat <- summarise(beh.all, acc = round(mean(resp.corr),2))
  beh.stat.rt <- summarise(beh.all[beh.all$resp.corr == 1, ], mean.rt = mean(resp.rt), sd.rt = sd(resp.rt), count = n())
  beh.stat <- data.frame(beh.stat, beh.stat.rt[, c('mean.rt', 'sd.rt', 'count')])
  beh.stat$se.rt <- beh.stat$sd.rt / sqrt(beh.stat$count)
  
  # count plot
  print('Plot pooled rating count...')
  p <- ggplot(data.rating.all, aes(x = factor(rating, levels = -5:5), fill = cond)) 
  p <- p + geom_bar(position = 'dodge')
  p <- p + scale_x_discrete(drop = FALSE)
  p <- p + labs(title = paste0('N=', length(subs)), x = 'Rating', y = 'Count')
  p <- p + scale_fill_manual('Block', labels = toupper(levels(data.rating$cond)), 
                             values = as.character(cond.colors[which(cond.colors[,'cond']== levels(data.rating$cond)), 'color']))
  p <- p + theme(axis.text.x = element_text(color = rating.colors))
  
  # mean score plot
  print('Plot individual rating count...')
  p2 <- ggplot(data.mean, aes(x = factor(sub), y = mean.score, fill = cond))
  p2 <- p2 + geom_bar(stat = 'identity', position = 'dodge')
  p2 <- p2 + labs(title = 'Mean score', x = 'Participant', y = 'Rating') 
  p2 <- p2 + scale_fill_manual('Block', labels = toupper(levels(data.rating$cond)), 
                             values = as.character(cond.colors[which(cond.colors[,'cond']== levels(data.rating$cond)), 'color']))
  # rating trend plot
  print('Plot pooled rating time course...')
  p3 <- ggplot(data.rating.all, aes(x = block, y = rating, col = cond))
  p3 <- p3 + geom_point(alpha = 0.5)  # in case of overlapping
  p3 <- p3 + geom_smooth(method = 'lm', se = FALSE) 
  p3 <- p3 + scale_y_continuous(limits = c(-5,5), breaks = -5:5)
  p3 <- p3 + labs(title = paste0('N=', length(subs)), x = 'Block sequence', y = 'Rating', caption = as.expression(bquote('Number in brackets indicates' ~ italic(R)^2)))
  p3 <- p3 + scale_color_manual('Block', labels = paste0(toupper(levels(data.rating$cond)), '(',temp$r2, temp$sig, ')'), 
                              values = as.character(cond.colors[which(cond.colors[,'cond']== levels(data.rating$cond)), 'color']))
  p3 <- p3 + theme(axis.text.y = element_text(color = rating.colors))
  
  # beh plot
  print('Plot individual behavioral data...')
  posn.d <- position_dodge(0.8)
  p4 <- ggplot(beh.stat, aes(x = factor(sub), y = mean.rt, fill = target))
  p4 <- p4 + geom_bar(position = posn.d, stat = 'identity')
  p4 <- p4 + geom_errorbar(aes(ymin = mean.rt - se.rt, ymax = mean.rt + se.rt), width = 0.4, position = posn.d)
  p4 <- p4 + labs(title = paste0('N=', length(subs)), x = 'Participant', y = 'Mean RT (s) and Accuracy')
  p4 <- p4 + geom_point(aes(x = factor(sub), y = acc), position = posn.d, shape = 9, size = 1, alpha = 1, color = 'blue')
  
  beh.stat  <- group_by(beh.stat, target)
  beh.group <- summarise(beh.stat, mean.rt.group = mean(mean.rt), acc.group = mean(acc), sd.rt = sd(mean.rt), sd.acc = sd(acc), count = n())
  beh.group$se.rt <- beh.group$sd.rt / sqrt(beh.group$count)
  beh.group$se.acc <- beh.group$sd.acc / sqrt(beh.group$count)

  print('Plot group average behavioral data...')
  p5 <- ggplot (beh.group, aes(x = target, y = mean.rt.group, fill = target))
  p5 <- p5 + geom_bar(position = posn.d, stat = 'identity')
  if (length(subs) > 1) {
    p5 <- p5 + geom_errorbar(aes(ymin = mean.rt.group - se.rt, ymax = mean.rt.group + se.rt), width = 0.3)
  }
  p5 <- p5 + labs(title = paste0('N=', length(subs)), x = 'Condition', y = 'Mean RT (s) and Accuracy')
  p5 <- p5 + geom_point(aes(x = target, y = acc.group), position = posn.d, shape = 9, size = 3, alpha = 1, color = 'blue')
  if (length(subs) > 1) {
    p5 <- p5 + geom_errorbar(aes(ymin = acc.group - se.acc, ymax = acc.group + se.acc), width = 0.2, color = 'blue')
  }
    
  if (savePlotOn){
    print('Ouput graphs...')
    ggarrange(p, p2, p3, nrow = 3)
    ggsave(filename = paste0('Rating result', subchar, '.tiff'), 
           path = f_result, device = 'tiff', width = 15, height = 30, units = 'cm')
    ggarrange(p4, p5, nrow = 2)
    ggsave(filename = paste0('Behavioral result', subchar, '.tiff'),
           path = f_result, device = 'tiff', width = 25, height = 20, units = 'cm')
  }
  
  if (individualOn) {
    print('Output graphs...')
    parrange <- get.plot.arrangement(length(subs))
    
    ggarrange(plotlist = plist, nrow = parrange[1], ncol = parrange[2])
    if (savePlotOn) {
      ggsave(filename = paste0('Rating result_indiv', subchar, '.tiff'), path = f_result, device = 'tiff', width = 40, height = 20, units = 'cm')
    }
    
    ggarrange(plotlist = glist, nrow = parrange[1], ncol = parrange[2])
    if (savePlotOn) {
      ggsave(filename = paste0('Rating trend_indiv', subchar, '.tiff'), path = f_result, device = 'tiff', width = 40, height = 20, units = 'cm')
    }
    
    ggarrange(plotlist = rlist, nrow = parrange[1], ncol = parrange[2])
    if (savePlotOn) {
      ggsave(filename = paste0('Beh_RT_indiv', subchar, '.tiff'), path = f_result, device = 'tiff', width = 40, height = 20, units = 'cm')
    }
    
    ggarrange(plotlist = alist, nrow = parrange[1], ncol = parrange[2])
    if (savePlotOn) {
      ggsave(filename = paste0('Beh_ACC_indiv', subchar, '.tiff'), path = f_result, device = 'tiff', width = 40, height = 20, units = 'cm')
    }
  }

}



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



get.data.sart <- function(sub, ureventOn = TRUE) {
  
  # path 
  f_subinfo <- file.path(f_main, 'subinfo.rdata')
  
  # pars
  cols2analyze <- c('number', 'type', 'trigger', 'blocks.thisN', 'trials.thisN', 
                    'resp.keys', 'resp.corr', 'resp.rt', 
                    'rating.response', 'rating.rt')
  
  data.raw <- read.csv(file.path(f_raw, paste0(sub, '_sart.csv')), header = TRUE, stringsAsFactors = FALSE)
  
  # register sub info
  register.subinfo(sub, data.raw)
  
  
  # filter out block_prac
  rmRowId   <- which(is.na(data.raw$rating.response) & is.na(data.raw$resp.corr))
  data <- data.raw[-rmRowId,]
  
  # gather columns
  data <- data[,cols2analyze]
  
  # add dis2pr
  data <- add.dis2pr(data)
  
  # add urevent
  if (ureventOn) {data <- add.urevent(data, sub)}
  
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



get.data.vs <- function(sub, ureventOn = TRUE) {
  

  # pars
  cols2analyze <- c('block', 'trigger', 'target', 'blocks.thisRepN', 'blocks.thisN', 'trials.thisN', 
                    'nTri', 'nSqr', 'nPen', 'nHex', 'resp.keys', 'resp.corr', 'resp.rt', 
                    'rating.response', 'rating.rt')
  
  data.raw <- read.csv(file.path(f_raw, paste0(sub, '_vs.csv')), header = TRUE, stringsAsFactors = FALSE)
  
  # register sub info
  register.subinfo(sub, data.raw)
  
  
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
  if (ureventOn) {
    data <- add.urevent(data, sub)
  }
  
  return(data)
}



add.urevent <- function(data, sub) {

  data$urevent <- NaN
  triggers <- unique(data$trigger)
  triggers <- triggers[!is.na(triggers)]
  
  temp <- readMat(file.path('urevent_matfile', paste0(sub, '.mat')))
  eeg.triggers <- temp[[1]]
  urevents <- eeg.triggers %in% triggers
  
  # consitency
  # in case the eeg recording was started after the behavior recording
  len2cut <- length(na.omit(data$trigger)) - length(eeg.triggers[urevents])
  if (len2cut == 0) {
    if (any(na.omit(data$trigger) != eeg.triggers[urevents])) {
      stop('Inconsidency between beh & eeg triggers!')
    }
  } else if (len2cut > 0) {
    if (any(na.omit(data[-c(1:len2cut),]$trigger) != eeg.triggers[urevents])) {
      stop('Inconsidency between beh & eeg triggers!')
    }
  } else {
    stop('Inconsidency between beh & eeg triggers!')
  }
  
  rows <- data$trigger %in% triggers
  if (len2cut > 0) {rows[1:len2cut] <- FALSE}
  data[rows,]$urevent <- which(urevents)
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



register.subinfo <- function(sub, data.raw) {
  
  f_subinfo <- file.path(f_main, 'subinfo.rdata')
  load(f_subinfo)
  newVec <- data.raw[1, c('participant', 'gender', 'age')]
  
  # check
  if (newVec[1] != sub) {
    warning('Sub id does not match the file record. Use the speficified sub id.')
    newVec[1] <- sub
  }
  
  if (sum(subinfo$sub == sub) == 0){ # data not exist, create one
    print(paste0('Register new data ', sub))
    subinfo[nrow(subinfo)+1, c('sub', 'gender', 'age')] <- newVec
  } else {# data exist, check
    print(paste0('Data ', sub, ' exist. Check consistency ...'))
    oldVec <- subinfo[subinfo$sub == sub,c('sub', 'gender', 'age')] 
    if (!all.equal(newVec, oldVec, check.names = FALSE, check.attributes = FALSE)){
      print('Subject information differs between csv files; keep the old vector!')
    } else {
      print('Consistent.')
    }
  }
  
  # update
  save(file = f_subinfo, subinfo)
}




get.plot.arrangement <- function(nPlot) {
  if (nPlot <= 2) {return(c(1, nPlot))}
  else if (nPlot <= 4) {return(c(2, 2))}
  else if (nPlot <= 9) {return(c(ceiling(nPlot/3), 3))}
  else if (nPlot <= 16) {return(c(ceiling(nPlot/4), 4))}
  else if (nPlot <= 25) {return(c(ceiling(nPlot/5), 5))}
  else if (nPlot <= 36) {return(c(ceiling(nPlot/6), 6))}
  else if (nPlot <= 49) {return(c(ceiling(nPlot/7), 7))}
}

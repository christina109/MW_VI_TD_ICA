
colors <- c('#87CEEB', '#FF6347')
scolors <- c(colors, "#00BFFF", "#FF4500")
names(scolors) <- c('ot', 'mw', 'eoa', 'ioa')


measures <- c('fAlpha', 'pALpha', 'loAlpha', 'roAlpha', 'fTheta', 'pTheta', 'loTheta', 'roTheta', 
              'loP1', 'roP1', 'loN1', 'roN1', 'P3', 
              'pfCohAlpha', 'ploCohAlpha', 'proCohAlpha', 'loroCohAlpha', 'lofCohAlpha', 'rofCohAlpha',
              'pfCohTheta', 'ploCohTheta', 'proCohTheta', 'loroCohTheta', 'lofCohTheta', 'rofCohTheta',
              paste0('alphaPC', 1:32), paste0('alphaOverallPC', 1:32),
              paste0('alphaOPC', 1:32, 'epochs'), paste0('alphaEOPC', 1:32, 'tNorm'),
              paste0('alphaEOICp', 1:32), paste0('alphaEOICo', 1:32))
feats    <- c('alpha', 'alpha', 'alpha', 'alpha', 'theta', 'theta', 'theta', 'theta',
              'p1', 'p1', 'n1', 'n1', 'p3',
              'alpha', 'alpha', 'alpha', 'alpha', 'alpha', 'alpha',
              'theta', 'theta', 'theta', 'theta', 'theta', 'theta',
              rep('alpha.pca', 32), rep('alpha.pca', 32), 
              rep('alpha.pca', 32), rep('alpha.pca', 32),
              rep('alpha.ica.conds.raw', 32), rep('alpha.ica.oscillation', 32))
folders  <- c('alpha_fz', 'alpha_pz', 'alpha_o1', 'alpha_o2', 
              'theta_fz', 'theta_pz', 'theta_o1', 'theta_o2', 
              'p1_o1', 'p1_o2', 'n1_o1', 'n1_o2', 'p3_pz',
              'ispc_alpha_fz_pz', 'ispc_alpha_pz_o1', 'ispc_alpha_pz_o2', 
              'ispc_alpha_o1_o2', 'ispc_alpha_fz_o1', 'ispc_alpha_fz_o2',
              'ispc_theta_fz_pz', 'ispc_theta_pz_o1', 'ispc_theta_pz_o2', 
              'ispc_theta_o1_o2', 'ispc_theta_fz_o1', 'ispc_theta_fz_o2',
              paste0('alpha_pca_', 1:32), paste0('alpha_pca_overall_', 1:32),
              paste0('alpha_pca_overall_', 1:32,'_epochs'), 
              paste0('alpha_pca_overall_', 1:32,'_epochs_tNorm'),
              paste0('alpha_ica_conds_raw_overall_', 1:32, '_epochs'),
              paste0('alpha_ica_oscillation_overall_', 1:32, '_epochs'))
measureNames <- c('alpha Fz', 'alpha Pz', 'alpha O1', 'alpha O2', 
                  'theta Fz', 'theta Pz', 'theta O1', 'theta O2', 
                  'P1 O1', 'P1 O2', 'N1 O1','N1 O2', 'P3 Pz',
                  'alpha Fz-Pz', 'alpha Pz-O1', 'alpha Pz-O2',
                  'alpha O1-O2',  'alpha Fz-O1', 'alpha Fz-O2',
                  'theta Fz-Pz', 'theta Pz-O1', 'theta Pz-O2',
                  'theta O1-O2',  'theta Fz-O1', 'theta Fz-O2',
                  paste0('alpha PC', 1:32), paste0('alpha overall PC', 1:32),
                  paste0('alpha overall PC', 1:32),  paste0('alpha overall PC', 1:32),
                  paste0('alpha overall IC', 1:32), paste0('alpha IC', 1:32))
measureNcol  <- c(rep(3, 4+4), rep(4, 5), rep(3, 6+6), rep(310, 32+32), 
                  rep(3,32), rep(3,32), rep(3, 32), rep(13, 32))
measureTypes <- c('power', 'power', 'power', 'power',
                  'power', 'power', 'power', 'power',
                  'ERP', 'ERP', 'ERP', 'ERP', 'ERP',
                  'ISPC', 'ISPC', 'ISPC', 'ISPC', 'ISPC', 'ISPC', 
                  'ISPC', 'ISPC', 'ISPC', 'ISPC', 'ISPC', 'ISPC',
                  rep('PCA', 32), rep('overall PCA', 32),
                  rep('overall PCA', 32), rep('overall PCA', 32),
                  rep('overall ICA', 32), rep('overall ICA', 32))
measureColors <- c("#FF0000", '#0000ff', '#00ff00', '#e1e11a',
                   "#FF0000", '#0000ff', '#00ff00', '#e1e11a',
                   '#00ff00', '#e1e11a', '#00ff00', '#e1e11a', '#0000ff',
                   '#7f007f', "#007f7f", "#7f7f7f", "#007f00", "#7f7f00", "#ff7f00", 
                   '#7f007f', "#007f7f", "#7f7f7f", "#007f00", "#7f7f00", "#ff7f00",
                   rainbow(32), rainbow(32),
                   rainbow(32), rainbow(32),
                   rep(rainbow(8),4), rep('black', 32)) 
measureShapes <- c(rep(21,8), rep(22,5), rep(24,12), rep(23,32), rep(23,32), 
                   rep(23,32), rep(23,32), rep(23,32), rep(23, 32))
# measureShaps <- c(rep(1,4), rep(2,4), rep(3,5), rep(1,6), rep(2,6))

hIncreaseOn <- as.logical(c(0, rep(1, 3), rep(0, 4 + 5 + 6 * 2 + 32 + 32 + 32 +32 + 32 + 32)))

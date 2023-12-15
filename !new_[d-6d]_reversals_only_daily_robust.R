# ==================================================================================================
# [6d] COMPUTE AND ANALYZE REVERSAL ESTMATION: robustness
#  - RUN TIME: 
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer); library(lfe); library(plm); library(extrafont)
library(sandwich); library(lmtest); library(car); library(plotrix)
library(Hmisc); library(estimatr); library(readxl)

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.55.0/bin/gswin64c.exe")

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')
source('./codes/ratfunc.R')
source('./codes/corstars.R')
source('./codes/starfunc.R')

# load data
load('./data_input/subtrace_with_smoothper.RData')
load('./data_input/vlm_robust.RData')
remove(submex_init)
gc()

# --------------------------------------------------------------------------------------------------
# (2) COMPUTE REVERSALS with alternative approaches
# --------------------------------------------------------------------------------------------------

cdeb <- submergent[, .(cusip_id, BOND_TYPE)][subtrace[, .N, cusip_id],, on = 'cusip_id']
cdeb <- cdeb[, c(1,2)]
subtrace <- cdeb[subtrace,, on = 'cusip_id']
#subtrace <- delete(subtrace, subtrace$BOND_TYPE == 'CDEB')
remove(cdeb)
#subtrace[, BOND_TYPE := NULL]
gc()

# compute market return
setkey(subtrace, cusip_id, trd_exctn_dt)
subtrace[, lsize := shift(size_hist), cusip_id]
subtrace <- delete(subtrace, keep.idxs = !is.na(subtrace$lsize))
subtrace[, mkt := weighted.mean(ret, lsize), trd_exctn_dt]
gc()

# [USUAL] define function to compute stuff
revfunc = function(r, v, iv) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 3))
    
  } else {
    
    m = lm('rprime ~ r + rv + riv', data = df)
    mod = summary(m)
    as.list(c(coef(mod)[2:4,1]))
    
  }
}

# [USUAL with half obs] define function to compute stuff
revfunc_fewerobs = function(r, v, iv) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < 0.25*Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 3))
    
  } else {
    
    m = lm('rprime ~ r + rv + riv', data = df)
    mod = summary(m)
    as.list(c(coef(mod)[2:4,1]))
    
  }
}

# [USUAL] define function to variance
revfunc_vars = function(r, v, iv) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 3))
    
  } else {
    
    m = lm('rprime ~ r + rv + riv', data = df)
    mod = summary(m)
    as.list(c(1/(coef(mod)[2:4,2])^2))

  }
}

# [PLUS MARKET] define function to compute stuff
revfunc_mkt = function(r, v, iv, mkt) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv, mkt = mkt)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 3))
    
  } else {
    
    mod = summary(lm('rprime ~ r + rv + riv + mkt', data = df))
    as.list(c(coef(mod)[2:4,1]))
    
  }
}

# [PLUS VOLUME] define function to compute stuff
revfunc_vlm = function(r, v, iv) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 3))
    
  } else {
    
    mod = summary(lm('rprime ~ r + rv + riv + v + iv', data = df))
    as.list(c(coef(mod)[2:4,1]))
    
  }
}

# [PLUS CUM VOLUME] define function to compute stuff
revfunc_cumvlm = function(r, v, iv, cumiv) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv, rcumiv = r*cumiv)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 3))
    
  } else {
    
    mod = summary(lm('rprime ~ r + rv + riv + rcumiv', data = df))
    as.list(c(coef(mod)[2:4,1]))
    
  }
}


# USUAL
subtrace[, c('C1_base', 'C2_base', 'C3_base') := revfunc(r = ret,
                                                         v = ctrnvr_alt_3,
                                                         iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# PRICE RETURN
subtrace[, c('C1_pret', 'C2_pret', 'C3_pret') := revfunc(r = plogret,
                                                         v = ctrnvr_alt_3,
                                                         iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()


# USUAL + log volumes
subtrace[, c('C1_logvlm', 'C2_logvlm', 'C3_logvlm') := revfunc(r = ret,
                                                               v = ctrnvr_alt_5,
                                                               iv = itrnvr_alt_5),
         .(cusip_id, ind)]
gc()

# USUAL + round-trip turnover
rt_ind <- subtrace[, .(.N, sum(rttrnvr > min(rttrnvr))), .(cusip_id, ind)][V2 > 10 & V2/N > 0.2]
rt_ind[, rt_ind := 1]
rt_ind <- rt_ind[, .(cusip_id, ind, rt_ind)]
subtrace <- rt_ind[subtrace,, on = c('cusip_id', 'ind')]
remove(rt_ind)
gc()
setkey(subtrace, cusip_id, trd_exctn_dt)

subtrace[rt_ind == 1, c('C1_rtvlm', 'C2_rtvlm', 'C3_rtvlm') := revfunc(r = ret,
                                                                       v = rttrnvr,
                                                                       iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# USUAL + VW mid price
subtrace[, c('C1_midp', 'C2_midp', 'C3_midp') := revfunc_fewerobs(r = ret_mid,
                                                                  v = ctrnvr_alt_3,
                                                                  iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# USUAL + large trades only
subtrace[, c('C1_bigp', 'C2_bigp', 'C3_bigp') := revfunc_fewerobs(r = ret_big,
                                                         v = ctrnvr_alt_3,
                                                         iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# USUAL + large trades above 100k only
subtrace[, c('C1_bigaltp', 'C2_bigaltp', 'C3_bigaltp') := revfunc_fewerobs(r = ret_big_alt,
                                                                     v = ctrnvr_alt_3,
                                                                     iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# USUAL + retail notes excluded
subtrace[!BOND_TYPE == 'RNT', c('C1_nornt', 'C2_nornt', 'C3_nornt') := revfunc(r = ret,
                                                                               v = ctrnvr_alt_3,
                                                                               iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# MKT 1st stage
subtrace[, c('C1_mkt', 'C2_mkt', 'C3_mkt') := revfunc_mkt(r = ret,
                                                          v = ctrnvr_alt_3,
                                                          iv = itrnvr_alt_3,
                                                          mkt = mkt),
         .(cusip_id, ind)]
gc()

# VLM in the 1st stage
subtrace[, c('C1_vlm', 'C2_vlm', 'C3_vlm') := revfunc_vlm(r = ret,
                                                          v = ctrnvr_alt_3,
                                                          iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# CUM VLM in the 1st stage
subtrace[, c('C1_cumvlm', 'C2_cumvlm', 'C3_cumvlm') := revfunc_cumvlm(r = ret,
                                                          v = ctrnvr_alt_3,
                                                          iv = itrnvr_alt_3,
                                                          cumiv = cumitrnvr_alt_4),
         .(cusip_id, ind)]
gc()

# DEALER VLM in the 1st stage
subtrace[, c('C1_dvlm', 'C2_dvlm', 'C3_dvlm') := revfunc_cumvlm(r = ret,
                                                                      v = ctrnvr_alt_3,
                                                                      iv = itrnvr_alt_3,
                                                                      cumiv = dtrnvr),
         .(cusip_id, ind)]
gc()

# vars of the base regressions
subtrace[, c('C1_var', 'C2_var', 'C3_var') := revfunc_vars(r = ret,
                                                           v = ctrnvr_alt_3,
                                                           iv = itrnvr_alt_3),
         .(cusip_id, ind)]
gc()

# --------------------------------------------------------------------------------------------------
# (3) CONSTRUCT CROSS-SECTION
# --------------------------------------------------------------------------------------------------

# name columns with coefficients
coefcols <- unlist(lapply(as.list(c('C1', 'C2', 'C3')), function(x) {
  paste0(x, c('_base', '_pret', '_logvlm', '_rtvlm', '_midp', '_bigp', '_bigaltp', '_nornt', '_mkt', '_vlm',
              '_cumvlm', '_dvlm', '_var'))
}))
  
# make a small summary
res = subtrace[, lapply(.SD, function(x) {x[1]}), .(cusip_id, issuer, ind, Nind),
               .SDcols = coefcols]

# filter
res = res[!is.na(C1_base) & !is.na(C2_base)& !is.na(C3_base)]

res = res[C1_base > quantile(C1_base, 0.01) & C1_base < quantile(C1_base, 0.99) &
          C2_base > quantile(C2_base, 0.01) & C2_base < quantile(C2_base, 0.99) &
          C3_base > quantile(C3_base, 0.01) & C3_base < quantile(C3_base, 0.99)]

# calculate average values of bond and issuer characteristics per active period
varcols <- c('size_hist', 'rating_num', 'stock_ba', 'SE', 'VE', 'dealers', 'nfunds', 'avg_ba')
toadd = subtrace[, lapply(.SD, function(z) {mean(z, na.rm =T)}), .(cusip_id, ind),
                 .SDcols = varcols]
res = toadd[res,,on = c('cusip_id', 'ind')]

# calculate initial values of bond and issuer characteristics per active period
setkey(subtrace, cusip_id, trd_exctn_dt)
toadd = subtrace[, lapply(.SD, function(z) {z[!is.na(z)][1]}), .(cusip_id, ind),
                 .SDcols = varcols]
names(toadd)[-c(1:2)] <- paste0(names(toadd)[-c(1:2)], '_init')
res = toadd[res,,on = c('cusip_id', 'ind')]

# create a cross-section of bonds by averaging across different active periods for each bond
# weighting with number of trading days in each active period
csres = res[, lapply(.SD, function(x) {weighted.mean(x, w = Nind, na.rm = T)}), .(cusip_id, issuer),
            .SDcols = c(coefcols, varcols)]

# add initial values from the very first active period
setkey(res, cusip_id, ind)
toadd <- res[, lapply(.SD, function(x) {x[1]}), .(cusip_id, issuer), .SDcols = c(varcols)]
names(toadd)[-c(1:2)] <- paste0(names(toadd)[-c(1:2)], '_init')
csres <- toadd[csres,, on = c('cusip_id', 'issuer')]

# add most frequent rating per bond
toadd = subtrace_prefilt[, .(rating_num[1], .N), c('cusip_id', 'ind')]
names(toadd)[3] = 'rating'
toadd = toadd[, rating[which.max(N)], cusip_id]
names(toadd)[2] = 'rating'
csres = toadd[csres,,on = 'cusip_id']

# add initial rating per bond
toadd = subtrace_prefilt[, .(rating_num[1], .N), c('cusip_id', 'ind')]
names(toadd)[3] = 'rating'
setkey(toadd, cusip_id, ind)
toadd = toadd[, rating[1], cusip_id]
names(toadd)[2] = 'rating_init'
csres = toadd[csres,,on = 'cusip_id']

# add average number of trading days
toadd = res[, mean(Nind, na.rm = T), cusip_id]
names(toadd)[2] = 'Nind'
csres = toadd[csres,,on = 'cusip_id']

csres[, rating_num := rating_num/100]
csres[, HY := as.numeric(rating_num*100 > 10)]
csres[, hynfunds := as.numeric(rating_num*100 > 10)*nfunds]
csres[, size_hist := size_hist/1000]

csres[, rating_num_init := rating_num_init/100]
csres[, HY_init := as.numeric(rating_num_init*100 > 10)]
csres[, hynfunds_init := as.numeric(rating_num_init*100 > 10)*nfunds_init]
csres[, size_hist_init := size_hist_init/1000]

# negative values 
csres[, avg_ba_neg := -avg_ba]
csres[, stock_ba_neg := -stock_ba]
csres[, nfunds_neg := -nfunds]
csres[, dealers_neg := -dealers]
csres[, VE_neg := -VE]
csres[, size_hist_neg := -size_hist]

csres[, avg_ba_neg_init := -avg_ba_init]
csres[, stock_ba_neg_init := -stock_ba_init]
csres[, nfunds_neg_init := -nfunds_init]
csres[, dealers_neg_init := -dealers_init]
csres[, VE_neg_init := -VE_init]
csres[, size_hist_neg_init := -size_hist_init]

# winsorize coefs for round-trip
csres[, C1_rtvlm := truncate(C1_rtvlm, 0.01, 0.99)]
csres[, C2_rtvlm := truncate(C2_rtvlm, 0.01, 0.99)]
csres[, C3_rtvlm := truncate(C3_rtvlm, 0.01, 0.99)]

# create principal components of info asymmetry

pc_func = function(csdata, vrbls) {
  
  # calculate principal components
  fml = formula(paste('~', paste(vrbls, collapse = '+')))
  pc = prcomp(formula = fml, data = csdata, center = T, scale = T, na.action = na.exclude)
  
  # return PC1 if all loadings are positive
  pc1 = pc$rotation[,1]
  flag <- sum(pc1 > 0) == length(pc1)
  if (sum(pc1 < 0) == length(pc1)) {flag = T; pc$x[,1] <- (-1)*pc$x[,1]} 
  if (flag) { A = pc$x[,1] } else { A = rep(NA_real_, nrow(csdata)) }
  return(as.numeric(A))
  
}

pc_cols_2 = c('avg_ba', 'nfunds_neg', 'dealers_neg', 'size_hist_neg')
pc_cols_2_init = paste0(pc_cols_2, '_init')

pc_cols_3 = c('nfunds_neg', 'dealers_neg', 'size_hist_neg')
pc_cols_3_init = paste0(pc_cols_3, '_init')

pc_cols = list(pc_cols_2, pc_cols_2_init, pc_cols_3, pc_cols_3_init)
pcs = lapply(pc_cols, function(x) {pc_func(csdata = csres, vrbls = x)})
names(pcs) = c('pc2', 'pc2_init', 'pc3', 'pc3_init')
pcs = as.data.table(pcs)

csres = cbind(csres, pcs)

# --------------------------------------------------------------------------------------------------
# (4) ESTIMATE 2nd stage (with PC2)
# --------------------------------------------------------------------------------------------------

# estimate with non-initial values and unmodified 2ns stage
estfunc = function(coef) {
  
  fmla <- formula(paste(coef, '~ pc2 | rating | 0 | rating'))
  m <- felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe')
  return(m)
  
}

# run estimations 
est = lapply(as.list(coefcols[!grepl('_var', coefcols)]), estfunc)
names(est) = coefcols[!grepl('_var', coefcols)]

# estimate with initial values and unmodified 2nd stage
estfunc_init = function(coef) {
  
  fmla <- formula(paste(coef, '~ pc2_init | rating_init | 0 | rating_init'))
  m <- felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe')
  return(m)
  
}

toadd <- lapply(as.list(c('C1_base', 'C2_base', 'C3_base')), estfunc_init)
names(toadd) <- paste0(c('C1_base', 'C2_base', 'C3_base'), '_init')
est <- c(est, toadd)

# estimate with non-initial values and modified 2nd stage
estfunc_wls = function(coef) {
  
  vec <- unlist(csres[, paste0(substr(coef, 1, 2), '_var'), with = F])
  fmla <- formula(paste(coef, '~ pc2 | rating | 0 | rating'))
  m <- felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe', weights = vec)
  return(m)
  
}

toadd <- lapply(as.list(c('C1_base', 'C2_base', 'C3_base')), estfunc_wls)
names(toadd) <- paste0(c('C1_base', 'C2_base', 'C3_base'), '_wls')
est <- c(est, toadd)

# PRINT OUT

estlist = as.list(c('_pret', '_logvlm', '_rtvlm', '_midp', '_bigp', '_bigaltp', '_nornt', '_mkt', '_vlm', '_cumvlm',
                    '_dvlm', '_base_init', '_base_wls'))

tab = lapply(estlist, function(i) {
  
  t(data.table(
    starfunc(x = as.numeric(sapply(est[names(est)[grepl(i, names(est))]],
                                   function(x) {summary(x)$coefficients[1]})),
             p = as.numeric(sapply(est[names(est)[grepl(i, names(est))]],
                                   function(x) {summary(x)$coefficients[4]})),
             th = 3),
    sefunc(x = as.numeric(sapply(est[names(est)[grepl(i, names(est))]],
                                 function(x) {summary(x)$coefficients[2]})),
           th = 3)
  ))
  
})

tab = do.call(rbind, tab)

# add volume-robust from the baseline file
tab = rbind(tab,
            starfunc(x = vlm_robust[1,1:3],
                     p = vlm_robust[3,1:3],
                     th = 3),
            sefunc(x = vlm_robust[2,1:3],
                   th = 3))

# call rows
tab <- cbind(c('Log-clean-price return', '~' ,
               'Log-volumes', '~' ,
               '1-hour roundtrip volumes', '~' ,
               'Avg. of VW buy and sell prices', '~',
               'Trades less than \\$10k excluded', '~',
               'Trades less than \\$100k excluded', '~',
               'Retail notes excluded', '~',
               'Market return added', '~',
               'Volumes added linearly', '~',
               'Lagged inventory added', '~',
               'Inter-dealer turnover added', '~',
               'PCs extracted from initial obs.', '~',
               'Weighted observations', '~',
               'Vlm. correlation controls', '~'),
             tab)
colnames(tab) = c('~', '$\\hat{\\beta}_1$', '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$')

addtorow = list()
addtorow$pos = list(0, 14, 22)
addtorow$command = c("\\hline \\multicolumn{4}{c}{A. Different inputs in the 1st stage} \\\\\n",
                     "\\hline \\multicolumn{4}{c}{B. Different models in the 1st stage} \\\\\n",
                     "\\hline \\multicolumn{4}{c}{C. Different 2nd stage} \\\\\n")

print(xtable(tab, auto = F,
             align = c("l|", "l|", rep('c', ncol(tab)-1)), type = "latex"),
      file = '../Apps/ShareLaTeX/Bond reversals/robustness_pc2.tex', floating = F, include.colnames = T,
      include.rownames = F, size = 'footnotesize',
      sanitize.text.function=function(x){x},
      hline.after = c(nrow(tab)),
      na.print = "", add.to.row = addtorow)

# --------------------------------------------------------------------------------------------------
# (4) ESTIMATE 2nd stage (with PC3)
# --------------------------------------------------------------------------------------------------

# estimate with non-initial values and unmodified 2ns stage
estfunc = function(coef) {
  
  fmla <- formula(paste(coef, '~ pc3 | rating | 0 | rating'))
  m <- felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe')
  return(m)
  
}

# run estimations 
est = lapply(as.list(coefcols[!grepl('_var', coefcols)]), estfunc)
names(est) = coefcols[!grepl('_var', coefcols)]

# estimate with initial values and unmodified 2nd stage
estfunc_init = function(coef) {
  
  fmla <- formula(paste(coef, '~ pc3_init | rating_init | 0 | rating_init'))
  m <- felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe')
  return(m)
  
}

toadd <- lapply(as.list(c('C1_base', 'C2_base', 'C3_base')), estfunc_init)
names(toadd) <- paste0(c('C1_base', 'C2_base', 'C3_base'), '_init')
est <- c(est, toadd)

# estimate with non-initial values and modified 2nd stage
estfunc_wls = function(coef) {
  
  vec <- unlist(csres[, paste0(substr(coef, 1, 2), '_var'), with = F])
  fmla <- formula(paste(coef, '~ pc3 | rating | 0 | rating'))
  m <- felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe', weights = vec)
  return(m)
  
}

toadd <- lapply(as.list(c('C1_base', 'C2_base', 'C3_base')), estfunc_wls)
names(toadd) <- paste0(c('C1_base', 'C2_base', 'C3_base'), '_wls')
est <- c(est, toadd)

# PRINT OUT

estlist = as.list(c('_pret', '_logvlm', '_rtvlm', '_midp', '_bigp', '_bigaltp', '_nornt', '_mkt', '_vlm', '_cumvlm',
                    '_dvlm', '_base_init', '_base_wls'))

tab = lapply(estlist, function(i) {
  
  t(data.table(
    starfunc(x = as.numeric(sapply(est[names(est)[grepl(i, names(est))]],
                                   function(x) {summary(x)$coefficients[1]})),
             p = as.numeric(sapply(est[names(est)[grepl(i, names(est))]],
                                   function(x) {summary(x)$coefficients[4]})),
             th = 3),
    sefunc(x = as.numeric(sapply(est[names(est)[grepl(i, names(est))]],
                                 function(x) {summary(x)$coefficients[2]})),
           th = 3)
  ))
  
})

tab = do.call(rbind, tab)

# add volume-robust from the baseline file
tab = rbind(tab,
            starfunc(x = vlm_robust[1,4:6],
                     p = vlm_robust[3,4:6],
                     th = 3),
            sefunc(x = vlm_robust[2,4:6],
                   th = 3))

# call rows
tab <- cbind(c('Log-clean-price return', '~' ,
               'Log-volumes', '~' ,
               '1-hour roundtrip volumes', '~' ,
               'Avg. of VW buy and sell prices', '~',
               'Trades less than \\$10k excluded', '~',
               'Trades less than \\$100k excluded', '~',
               'Retail notes excluded', '~',
               'Market return added', '~',
               'Volumes added linearly', '~',
               'Lagged inventory added', '~',
               'Inter-dealer turnover added', '~',
               'PCs extracted from initial obs.', '~',
               'Weighted observations', '~',
               'Vlm. correlation controls', '~'),
             tab)
colnames(tab) = c('~', '$\\hat{\\beta}_1$', '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$')

addtorow = list()
addtorow$pos = list(0, 14, 22)
addtorow$command = c("\\hline \\multicolumn{4}{c}{A. Different inputs in the 1st stage} \\\\\n",
                     "\\hline \\multicolumn{4}{c}{B. Different models in the 1st stage} \\\\\n",
                     "\\hline \\multicolumn{4}{c}{C. Different 2nd stage} \\\\\n")

print(xtable(tab, auto = F,
             align = c("l|", "l|", rep('c', ncol(tab)-1)), type = "latex"),
      file = '../Apps/ShareLaTeX/Bond reversals/robustness_pc3.tex', floating = F, include.colnames = T,
      include.rownames = F, size = 'footnotesize',
      sanitize.text.function=function(x){x},
      hline.after = c(nrow(tab)),
      na.print = "", add.to.row = addtorow)
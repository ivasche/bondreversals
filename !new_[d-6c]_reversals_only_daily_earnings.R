# ==================================================================================================
# [6c] COMPUTE AND ANALYZE REVERSAL ESTMATION: account for days to earnings announcement
#  - RUN TIME: 
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer); library(lfe); library(plm); library(extrafont)
library(sandwich); library(lmtest); library(car); library(plotrix)
library(Hmisc); library(estimatr)

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.55.0/bin/gswin64c.exe")

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')
source('./codes/ratfunc.R')
source('./codes/corstars.R')

# load data
load('./data_input/subtrace_with_smoothper.RData')

# --------------------------------------------------------------------------------------------------
# (2) COMPUTE REVERSALS
# --------------------------------------------------------------------------------------------------

# select only traded part with earnings announcements
subtrace_ea <- subtrace[!is.na(distea)]

# mark periods
subtrace_ea[, eaind := as.numeric(distea >= 1 & distea <= 10)]

# filter smooth periods with at least 10 days in eaind
subtrace_ea[, Nea := sum(eaind), .(cusip_id, ind)]
subtrace_ea <- subtrace_ea[Nea >= 10]

# make a variable for pos and neg ea
subtrace_ea[, eaposind := sum(unique(sign(suescore)), na.rm = T), .(cusip_id, ind)]

# define function to compute stuff
revfunc = function(r, v, iv, ivsign, ea) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv,
                  rvea = r*v*ea, rivea = r*iv*ea)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 12))
    
  } else {
    
    m = lm('rprime ~ r + rv + riv + rvea + rivea', data = df)
    mod = summary(m)
    A <- as.list(c(coef(mod)[2:6,1], coef(mod)[2:6,3]))
    return(A)
    
  }
}

# DAILY --------------------------------------------------------------------------------------------

# apply it
subtrace_ea[, c('C1', 'C2', 'C3', 'C4', 'C5',
                't1', 't2', 't3', 't4', 't5') := revfunc(r = ret, v = ctrnvr_alt_3, iv = itrnvr_alt_3, ea = eaind),
            .(cusip_id, ind)]
gc()

# make a small summary
res = subtrace_ea[, .(C1[1], C2[1], C3[1], C4[1], C5[1], t1[1], t2[1], t3[1], t4[1], t5[1]),
                  .(cusip_id, issuer, ind, Nind, eaposind)]
names(res)[6:ncol(res)] = c(paste0('C', 1:5), paste0('t', 1:5))

# filter
res = res[!is.na(C1) & !is.na(C2)& !is.na(C3) & !is.na(C4) & !is.na(C5)]

res = res[C1 > quantile(C1, 0.01) & C1 < quantile(C1, 0.99) &
            C2 > quantile(C2, 0.01) & C2 < quantile(C2, 0.99) &
            C3 > quantile(C3, 0.01) & C3 < quantile(C3, 0.99) & 
            C4 > quantile(C4, 0.01) & C4 < quantile(C4, 0.99) &
            C5 > quantile(C5, 0.01) & C5 < quantile(C5, 0.99)]

# means
round(res[, lapply(.SD, function(x) {coef(summary(lm(x ~ 1)))[c(1, 4)]}),,
          .SDcols = paste0('C', 1:5)], 3)

# medians
round(res[, lapply(.SD, median),, .SDcols = paste0('C', 1:5)], 3)

# print summary
tab = rbind(res[, .(mean(C1), median(C1), sum(C1>0), sum(C1<0),
                    sum((C1>0)*(t1>1.65)), sum((C1<0)*(t1 < -1.65)), length(C1))],
            res[, .(mean(C2), median(C2), sum(C2>0), sum(C2<0),
                    sum((C2>0)*(t2>1.65)), sum((C2<0)*(t2 < -1.65)), length(C2))],
            res[, .(mean(C3), median(C3), sum(C3>0), sum(C3<0),
                    sum((C3>0)*(t3>1.65)), sum((C3<0)*(t3 < -1.65)), length(C3))],
            res[, .(mean(C4), median(C4), sum(C4>0), sum(C4<0),
                    sum((C4>0)*(t4>1.65)), sum((C4<0)*(t4 < -1.65)), length(C4))],
            res[, .(mean(C5), median(C5), sum(C5>0), sum(C5<0),
                    sum((C5>0)*(t5>1.65)), sum((C5<0)*(t5 < -1.65)), length(C5))])
colnames(tab) = c('~Mean~', '~Med.~', 'No.$>$0~', 'No.$<$0~', 'No.$>$0*', 'No.$<$0*', 'No. Obs.')
rownames(tab) = c('$\\hat{\\beta}_1$','$\\hat{\\beta}_2$','$\\hat{\\beta}_3$',
                  '$\\hat{\\beta}_4$', '$\\hat{\\beta}_5$')

print(xtable(tab, align = c('l|', 'c', 'c|', 'c', 'c|', 'c', 'c|', 'c'),
             digits = c(0, 4, 4, rep(0, 5)), auto = F),
      file = '../Apps/ShareLaTeX/Bond reversals/restab_ea.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)))

# --------------------------------------------------------------------------------------------------
# (3) ANALYZE REVERSALS
# --------------------------------------------------------------------------------------------------

# calculate average values of bond and issuer characteristics per active period
toadd = subtrace_ea[, lapply(.SD, function(z) {mean(z, na.rm =T)}), .(cusip_id, ind),
                 .SDcols = c('size_hist', 'rating_num', 'stock_ba', 'SE', 'VE',
                             'dealers', 'nfunds', 'avg_ba')]
res = toadd[res,,on = c('cusip_id', 'ind')]

# create a cross-section of bonds by averaging across different active periods for each bond
# weighting with number of trading days in each active period
csres = res[, lapply(.SD, function(x) {weighted.mean(x, w = Nind, na.rm = T)}), .(cusip_id, issuer, eaposind),
            .SDcols = c('C1', 'C2', 'C3', 'C4', 'C5',
                        'size_hist', 'rating_num', 'stock_ba', 'SE', 'VE',
                        'dealers', 'nfunds', 'avg_ba')]

# add most frequent rating per bond
toadd = subtrace_prefilt[, .(rating_num[1], .N), c('cusip_id', 'ind')]
names(toadd)[3] = 'rating'
toadd = toadd[, rating[which.max(N)], cusip_id]
names(toadd)[2] = 'rating'
csres = toadd[csres,,on = 'cusip_id']

# add average number of trading days
toadd = res[, mean(Nind, na.rm = T), cusip_id]
names(toadd)[2] = 'Nind'
csres = toadd[csres,,on = 'cusip_id']

csres[, rating_num := rating_num/100]
csres[, HY := as.numeric(rating_num*100 > 10)]
csres[, hynfunds := as.numeric(rating_num*100 > 10)*nfunds]
csres[, size_hist := size_hist/1000]

# negative values 
csres[, avg_ba_neg := -avg_ba]
csres[, stock_ba_neg := -stock_ba]
csres[, nfunds_neg := -nfunds]
csres[, dealers_neg := -dealers]
csres[, VE_neg := -VE]
csres[, size_hist_neg := -size_hist]

# create principal components of info asymmetry ----------------------------------------------------

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

pc_cols_1 = c('avg_ba', 'nfunds_neg', 'dealers_neg', 'size_hist_neg', 'VE_neg', 'stock_ba')
pc_cols_2 = c('avg_ba', 'nfunds_neg', 'dealers_neg', 'size_hist_neg')
pc_cols_3 = c('nfunds_neg', 'dealers_neg', 'size_hist_neg')

pc_cols = list(pc_cols_1, pc_cols_2, pc_cols_3)
pcs = lapply(pc_cols, function(x) {pc_func(csdata = csres, vrbls = x)})
names(pcs) = paste0('pc', 1:length(pcs))
pcs = as.data.table(pcs)

csres = cbind(csres, pcs)

# run some summary stats
sumtabfunc = function(x) {
  
  tab = rbind(descfull(x[, C1]), descfull(x[, C2]), descfull(x[, C3]), descfull(x[, C4]),
              descfull(x[, C5]),
              descfull(x[, rating]),
              descfull(x[, avg_ba]),  
              descfull(x[, nfunds], dig = 1),
              descfull(x[, size_hist]),
              descfull(x[, dealers], dig = 1),
              descfull(x[, VE], dig = 1),
              descfull(x[, stock_ba]),
              descfull(x[, pc1]),
              descfull(x[, pc2]),
              descfull(x[, pc3])
  )
  
  rownames(tab) = c('$\\hat{\\beta}_1$', '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$',
                    '$\\hat{\\beta}_4$', '$\\hat{\\beta}_5$',
                    'Credit rating',
                    'Bond bid-ask, \\%',
                    'No. mutual fund owners',
                    'Issue size, bln \\$',
                    'No. dealers',
                    'Issuer size, bln \\$',
                    'Stock bid-ask, \\%',
                    'PC$_\\text{all}$',
                    'PC$_\\text{bond}$',
                    'PC$_\\text{bond-ex-ba}$'
  )
  tab
  
}

tab = sumtabfunc(csres)

print(xtable(tab, align = c('l|', rep('c', 2), 'c|', rep('c', ncol(tab)-4), '|c'),
             auto = F), file = '../Apps/ShareLaTeX/Bond reversals/sumstat_cs_D_ea.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, 5, 6, 12, nrow(tab)))

sdnms = c('avg_ba', 'avg_ba_neg', 'nfunds', 'nfunds_neg', 'size_hist', 'size_hist_neg',
          'dealers', 'dealers_neg', 'VE', 'VE_neg', 'stock_ba', 'stock_ba_neg',
          'rating_num')
csres[, c('avg_ba', 'avg_ba_neg', 'nfunds', 'nfunds_neg', 'size_hist', 'size_hist_neg',
          'dealers', 'dealers_neg', 'VE', 'VE_neg', 'stock_ba', 'stock_ba_neg',
          'rating_num') := lapply(.SD, function(x) {x/sd(x, na.rm = T)}),
      .SDcols = sdnms]

# run regressions of C1, C2, and C3 ----------------------------------------------------------------

estfunc = function(coef) {
  
  fmla = formula(paste(coef, '~ pc2 | rating | 0 | rating'))
  
  m = felm(fmla, data = as.data.frame(csres), keepX = T, cmethod='reghdfe')
  return(m)
  
}

# run estimations 
est_ea = lapply(as.list(paste0('C', 1:5)), estfunc)

# save joint
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc_ea.tex',
          est_ea,
          font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
          covariate.labels = c('PC$_\\text{bond}$'),
          float = F, intercept.top = T, intercept.bottom = F, digits = 3,
          dep.var.labels.include = T, model.numbers = F,
          dep.var.labels = c('$\\hat{\\beta}$\\textsubscript 1',
                             '$\\hat{\\beta}$\\textsubscript 2',
                             '$\\hat{\\beta}$\\textsubscript 3',
                             '$\\hat{\\beta}$\\textsubscript 4',
                             '$\\hat{\\beta}$\\textsubscript 5'),
          dep.var.caption = '',
          add.lines = list(c("Rating FE", rep('YES', length(est_ea)))),
          notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc_ea.tex')

# run regressions of C1, C2, and C3 splitting on positivity ------------------------------------

estfunc_posneg = function(coef) {
  
  m = vector('list', 2)
  
  fmla = formula(paste(coef, '~ pc2 | rating | 0 | rating'))
  
  m[[1]] = felm(fmla, data = as.data.frame(csres[eaposind==1]), keepX = T, cmethod='reghdfe')
  m[[2]] = felm(fmla, data = as.data.frame(csres[eaposind==-1]), keepX = T, cmethod='reghdfe')
  
  return(m)
  
}

# run estimations 
est_ea_posneg = lapply(as.list(paste0('C', 1:5)), estfunc_posneg)
est_ea_posneg = unlist(est_ea_posneg, recursive = F)

# save joint
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc_ea_posneg.tex',
            est_ea_posneg,
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = T, model.numbers = F,
            column.labels = rep(c('Pos', 'Neg'), 5),
            dep.var.labels = c('$\\hat{\\beta}_1$', '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$',
                               '$\\hat{\\beta}_4$', '$\\hat{\\beta}_5$'),
            dep.var.caption = '',
            add.lines = list(c("Rating FE", rep('YES', length(est_ea_posneg)))),
            notes.align = 'l')
)
mod <- sub("lcccccccccc", "lcc|cc|cc|cc|cc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc_ea_posneg.tex')
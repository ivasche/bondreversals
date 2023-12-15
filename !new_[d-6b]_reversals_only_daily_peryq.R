# ==================================================================================================
# [6v2_only_daily] COMPUTE AND ANALYZE REVERSAL ESTMATION
#  - RUN TIME: 
# ==================================================================================================

# load fonts
library(extrafont)
library(ggplot2)

# load data.table package
library(Hmisc)
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer); library(lfe); library(plm); 
library(sandwich); library(lmtest); library(car); library(plotrix)

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.55.0/bin/gswin64c.exe")

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')
source('./codes/ratfunc.R')
source('./codes/corstars.R')

# load data
load('./data_input/subtrace_with_smoothper_peryq.RData')

# --------------------------------------------------------------------------------------------------
# (2) COMPUTE REVERSALS
# --------------------------------------------------------------------------------------------------

# define function to compute stuff
revfunc = function(r, v, iv, ivsign) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, iv = iv, rv = r*v, riv = r*iv, ivsign = ivsign)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 12))
    
  } else {
    
    m = lm('rprime ~ r + rv + riv', data = df)
    mod = summary(m)
    as.list(c(coef(mod)[2:4,1],
              coef(mod)[2:4,3],
              coef(summary(lm(v[2:length(v)] ~ v[1:(length(v)-1)])))[2, c(1, 4)],
              coef(summary(lm(iv[2:length(iv)] ~ iv[1:(length(iv)-1)])))[2, c(1, 4)],
              c(cor.test(v, iv)$estimate, cor.test(v, iv)$p.value),
              c(cor.test(v, ivsign)$estimate, cor.test(v, ivsign)$p.value),
              as.numeric(coef(mod)[3, 1] - coef(mod)[4,1] > 0),
              linearHypothesis(m, 'rv = riv')$F[2]))
    
  }
}

# DAILY --------------------------------------------------------------------------------------------

# apply it
subtrace[, ind := NULL]
subtrace[, c('C1', 'C2', 'C3', 't1', 't2', 't3', 'corv', 'pcorv',
             'delta', 'pdelta', 'corviv', 'pcorviv', 'corvivsign', 'pcorvivsign',
             'coefdiff', 'pcoefdiff') := revfunc(r = ret, v = ctrnvr_alt_3, iv = itrnvr_alt_3,
                                                 ivsign = itrnvr_alt_4),
         .(cusip_id, yq)]
gc()

# make a small summary
res = subtrace[, .(C1[1], C2[1], C3[1], t1[1], t2[1], t3[1], corv[1], pcorv[1],
                   delta[1], pdelta[1], corviv[1], pcorviv[1], corvivsign[1], pcorvivsign[1],
                   coefdiff[1], pcoefdiff[1]), .(cusip_id, issuer, yq, Nind)]
names(res)[5:20] = c('C1', 'C2', 'C3', 't1', 't2', 't3',
                     'corv', 'pcorv', 'delta', 'pdelta', 'corviv', 'pcorviv',
                     'corvivsign', 'pcorvivsign', 'coefdiff', 'pcoefdiff')

# filter
res = res[!is.na(C1) & !is.na(C2)& !is.na(C3)]

res = res[C1 > quantile(C1, 0.01) & C1 < quantile(C1, 0.99) &
            C2 > quantile(C2, 0.01) & C2 < quantile(C2, 0.99) &
            C3 > quantile(C3, 0.01) & C3 < quantile(C3, 0.99)]

# means
round(res[, lapply(.SD, function(x) {coef(summary(lm(x ~ 1)))[c(1, 4)]}),,
          .SDcols = c('C1', 'C2', 'C3', 'corv', 'delta', 'corviv')], 3)

# medians
round(res[, lapply(.SD, median),, .SDcols = c('C1', 'C2', 'C3', 'corv', 'delta', 'corviv')], 3)

# print summary
tab = rbind(res[, .(mean(C1), median(C1), sum(C1>0), sum(C1<0),
                    sum((C1>0)*(t1>1.65)), sum((C1<0)*(t1 < -1.65)), length(C1))],
            res[, .(mean(C2), median(C2), sum(C2>0), sum(C2<0),
                    sum((C2>0)*(t2>1.65)), sum((C2<0)*(t2 < -1.65)), length(C2))],
            res[, .(mean(C3), median(C3), sum(C3>0), sum(C3<0),
                    sum((C3>0)*(t3>1.65)), sum((C3<0)*(t3 < -1.65)), length(C3))],
            res[, .(mean(C2-C3), median(C2-C3), sum(C2-C3>0), sum(C2-C3<0),
                    sum((C2-C3>0)*(pcoefdiff>0.10)), sum((C2-C3<0)*(pcoefdiff>0.10)), length(C3))])
colnames(tab) = c('~Mean~', '~Med.~', 'No.$>$0~', 'No.$<$0~', 'No.$>$0*', 'No.$<$0*', 'No. Obs.')
rownames(tab) = c('$\\hat{\\beta}_1$','$\\hat{\\beta}_2$','$\\hat{\\beta}_3$',
                  '$\\hat{\\beta}_2 - \\hat{\\beta}_3$')

print(xtable(tab, align = c('l|', 'c', 'c|', 'c', 'c|', 'c', 'c|', 'c'),
             digits = c(0, 4, 4, rep(0, 5)), auto = F),
      file = '../Apps/ShareLaTeX/Bond reversals/restab_peryq.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)))

# print volume correlations
tab = rbind(res[, .(mean(corviv), median(corviv), sum(corviv>0), sum(corviv<0),
                    sum((corviv>0)*(pcorviv<0.10)), sum((corviv<0)*(pcorviv < 0.10)), length(pcorviv))],
            res[, .(mean(corvivsign), median(corvivsign), sum(corvivsign>0), sum(corvivsign<0),
                    sum((corvivsign>0)*(pcorvivsign<0.10)), sum((corvivsign<0)*(pcorvivsign < 0.10)),
                    length(pcorvivsign))],
            res[, .(mean(corv), median(corv), sum(corv>0), sum(corv<0),
                    sum((corv>0)*(pcorv<0.10)), sum((corv<0)*(pcorv < 0.10)), length(corv))],
            res[, .(mean(delta), median(delta), sum(delta>0), sum(delta<0),
                    sum((delta>0)*(pdelta<0.10)), sum((delta<0)*(pdelta < 0.10)), length(delta))])
colnames(tab) = c('~Mean~', '~Med.~', 'No.$>$0~', 'No.$<$0~', 'No.$>$0*', 'No.$<$0*', 'No. Obs.')
rownames(tab) = c('$\\text{Corr}(V^{(c)}_t, \\lvert V^{(s)}_t \\rvert)$',
                  '$\\text{Corr}(V^{(c)}_t, V^{(s)}_t)$',
                  '$\\text{Corr}(V^{(c)}_{t}, V^{(c)}_{t-1})$',
                  '$\\text{Corr}(\\lvert V^{(s)}_{t} \\rvert, \\lvert V^{(s)}_{t-1} \\rvert)$')

print(xtable(tab, align = c('l|', 'c', 'c|', 'c', 'c|', 'c', 'c|', 'c'),
             digits = c(0, 3, 3, rep(0, 5)), auto = F),
      file = '../Apps/ShareLaTeX/Bond reversals/corvol_peryq.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)))

# --------------------------------------------------------------------------------------------------
# (3) ANALYZE REVERSALS
# --------------------------------------------------------------------------------------------------

# calculate average values of bond and issuer characteristics per active period
toadd = subtrace[, lapply(.SD, function(z) {mean(z, na.rm =T)}), .(cusip_id, yq),
                 .SDcols = c('size_hist', 'rating_num', 'stock_ba', 'SE', 'VE',
                             'dealers', 'nfunds', 'avg_ba')]
res = toadd[res,,on = c('cusip_id', 'yq')]

# average out so that there're one observation per bond per yq
csres = res

# add most frequent rating per bond per yq
toadd = subtrace_prefilt[, .(rating_num[1], .N), c('cusip_id', 'yq')]
names(toadd)[3] = 'rating'
toadd = toadd[, rating[which.max(N)], cusip_id]
names(toadd)[2] = 'rating'
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

# clean
remove(toadd, toadd_mean, tab)

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
  
  tab = rbind(descfull(x[, C1]), descfull(x[, C2]), descfull(x[, C3]),
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
             auto = F), file = '../Apps/ShareLaTeX/Bond reversals/sumstat_cs_yq.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, 3, 4, 10, nrow(tab)))

tab = corstars(csres[, .(avg_ba, nfunds, size_hist, dealers, VE, stock_ba)])
colnames(tab) =  c('Bond bid-ask', 'No. funds', 'Issue size', 'No. dealers', 'Issuer size')
rownames(tab) =  c('Bond bid-ask', 'No. funds', 'Issue size', 'No. dealers', 'Issuer size', 'Stock bid-ask')

print(xtable(tab[2:6,], align = c('l|', rep('c', 5)),
             auto = F), file = '../Apps/ShareLaTeX/Bond reversals/sumcorr_cs_yq.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)-1))

# standardize
sdnms = c('avg_ba', 'avg_ba_neg', 'nfunds', 'nfunds_neg', 'size_hist', 'size_hist_neg',
          'dealers', 'dealers_neg', 'VE', 'VE_neg', 'stock_ba', 'stock_ba_neg',
          'corv', 'corviv', 'rating_num')
csres[, c('avg_ba', 'avg_ba_neg', 'nfunds', 'nfunds_neg', 'size_hist', 'size_hist_neg',
          'dealers', 'dealers_neg', 'VE', 'VE_neg', 'stock_ba', 'stock_ba_neg',
          'corv', 'corviv', 'rating_num') := lapply(.SD, function(x) {x/sd(x, na.rm = T)}),
      .SDcols = sdnms]

# run regressions of C1, C2, and C3 ----------------------------------------------------------------

estfunc = function(coef) {
  
  fmla = vector('list', 14)
  
  fmla[[1]] = formula(paste(coef, '~ avg_ba | rating + yq | 0 | rating + yq'))
  fmla[[2]] = formula(paste(coef, '~ nfunds_neg | rating + yq | 0 | rating + yq'))
  fmla[[3]] = formula(paste(coef, '~ size_hist_neg | rating + yq | 0 | rating + yq'))
  fmla[[4]] = formula(paste(coef, '~ dealers_neg | rating + yq | 0 | rating + yq'))
  fmla[[5]] = formula(paste(coef, '~ VE_neg | rating + yq | 0 | rating + yq'))
  fmla[[6]] = formula(paste(coef, '~ stock_ba | rating + yq | 0 | rating + yq'))
  
  fmla[[7]] = formula(paste(coef, '~ avg_ba + nfunds_neg + size_hist_neg + dealers_neg | rating + yq | 0 | rating + yq'))
  fmla[[8]] = formula(paste(coef, '~ avg_ba + nfunds_neg + size_hist_neg + dealers_neg + VE_neg + stock_ba | rating + yq | 0 | rating + yq'))
  
  fmla[[9]] = formula(paste(coef, '~ pc1 | rating + yq | 0 | rating + yq'))
  fmla[[10]] = formula(paste(coef, '~ pc2 | rating + yq | 0 | rating + yq'))
  fmla[[11]] = formula(paste(coef, '~ pc3 | rating + yq | 0 | rating + yq'))
  
  fmla[[12]] = formula(paste(coef, '~ corv + corviv + pc1 | rating + yq | 0 | rating + yq'))
  fmla[[13]] = formula(paste(coef, '~ corv + corviv + pc2 | rating + yq | 0 | rating + yq'))
  fmla[[14]] = formula(paste(coef, '~ corv + corviv + pc3 | rating + yq | 0 | rating + yq'))
  
  mlist = lapply(fmla, function(x) {felm(x, data = as.data.frame(csres), keepX = T,
                                         cmethod='reghdfe')})
}

# run estimations 
est_c1   = estfunc(coef = 'C1')
est_c2 = estfunc(coef = 'C2')
est_c3 = estfunc(coef = 'C3')

# save output for c1
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_yq.tex', est_c1[1:8],
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('$\\text{Bond bid-ask}$',
                                 '$\\text{--No. funds}$',
                                 '$\\text{--Issue size}$',
                                 '$\\text{--No. dealers}$',
                                 '$\\text{--Issuer size}$',
                                 '$\\text{Stock bid-ask}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = F, model.numbers = T,
            dep.var.caption = '',
            add.lines = list(c("Rating FE", rep('YES', length(est_c1[1:8]))),
                             c("Time FE", rep('YES', length(est_c1[1:8])))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_yq.tex')
  
# save output for c2
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test2_yq.tex', est_c2[1:8],
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('$\\text{Bond bid-ask}$',
                                 '$\\text{--No. funds}$',
                                 '$\\text{--Issue size}$',
                                 '$\\text{--No. dealers}$',
                                 '$\\text{--Issuer size}$',
                                 '$\\text{Stock bid-ask}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = F, model.numbers = T,
            dep.var.caption = 'Dependent variable: $\\hat{\\beta}_2$',
            add.lines = list(c("Rating FE", rep('YES', length(est_c1[1:8]))),
                             c("Time FE", rep('YES', length(est_c1[1:8])))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test2_yq.tex')

# save output for c3
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test3_yq.tex', est_c3[1:8],
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('$\\text{Bond bid-ask}$',
                                 '$\\text{--No. funds}$',
                                 '$\\text{--Issue size}$',
                                 '$\\text{--No. dealers}$',
                                 '$\\text{--Issuer size}$',
                                 '$\\text{Stock bid-ask}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = F, model.numbers = T,
            dep.var.caption = 'Dependent variable: $\\hat{\\beta}_3$',
            add.lines = list(c("Rating FE", rep('YES', length(est_c1[1:8]))),
                             c("Time FE", rep('YES', length(est_c1[1:8])))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test3_yq.tex')

# save joint
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc_yq.tex',
            c(est_c1[10], est_c2[10], est_c3[10]),
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('\\hline PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = F, model.numbers = F,
            column.labels = c('$\\hat{\\beta}$\\textsubscript 1',
                              '$\\hat{\\beta}$\\textsubscript 2',
                              '$\\hat{\\beta}$\\textsubscript 3'),
            column.separate = NULL,          
            dep.var.caption = '',
            add.lines = list(c("Rating FE", rep('YES', 3)),
                             c("Time FE", rep('YES', 3))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc_yq.tex')


# RUN SEPARATELY PRE- AND POST-GFC

estfunc_prepost = function(coef, prepost) {
  
  fmla = formula(paste(coef, '~ pc2 | rating + yq | 0 | rating + yq'))

  if (prepost == 'pre') {
    csres_in = csres[yq <= as.yearqtr('2008 Q2')] 
  } else {
    csres_in = csres[yq >= as.yearqtr('2010 Q1')] 
  }
  
  m = felm(fmla, data = as.data.frame(csres_in), keepX = T, cmethod='reghdfe')
  return(m)
  
}

est_prepost = lapply(as.list(as.data.table(t(expand.grid(c('pre', 'post'), c('C1', 'C2', 'C3'))))),
                     function(i) {
                       estfunc_prepost(coef = i[2], prepost = i[1])
                     })

# save output
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc_yq_prepost.tex',
            est_prepost,
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = T, model.numbers = F,
            column.labels = rep(c('Pre-GFC', 'Post-GFC'), 3),
            dep.var.caption = '',
            dep.var.labels = c('$\\hat{\\beta}_1$', '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$'),
            add.lines = list(c("Rating FE", rep('YES', length(est_prepost)))),
            notes.align = 'l')
)
mod <- sub("lcccccc", "lcc|cc|cc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc_yq_prepost.tex')

# plot fixed effects

fedf = rbind(cbind(as.data.table(getfe(est_c1[[10]]))[fe == 'yq', .(idx, effect)], 'C1'),
             cbind(as.data.table(getfe(est_c2[[10]]))[fe == 'yq', .(idx, effect)], 'C2'),
             cbind(as.data.table(getfe(est_c3[[10]]))[fe == 'yq', .(idx, effect)], 'C3'))
names(fedf)[3] = 'coef'
fedf$coef <- as.factor(fedf$coef)
levels(fedf$coef) <- c(expression(hat(beta)[1]),
                       expression(hat(beta)[2]),
                       expression(hat(beta)[3]))

fedf[, idx := as.Date(as.yearmon(as.yearqtr(as.character(idx)))+3/12)-1]

p <- fedf %>%
  ggplot(mapping = aes(x = idx, y = effect))+
  geom_line(lwd = 2) + 
  #geom_point(size = 3, shape = 21, color = 'black', fill = 'white', stroke = 3) + 
  guides(color = "none", fill = "none") + 
  facet_wrap(~ as.factor(coef), nrow = 3, scales = 'free_y', labeller = label_parsed) +
  theme_bw(base_size = 26) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "lightgray"), 
        panel.grid.minor = element_blank(), axis.title = element_blank(), 
        text=element_text(family="CM Roman")) + 
  #geom_vline(xintercept = -5, linetype = "dotted", color = "red", lwd = 1.5) +
  #geom_vline(xintercept = 0, linetype = "dotted", color = "blue", lwd = 1.5) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")
p

ggsave(plot = p, filename = '../Apps/ShareLaTeX/Bond reversals/betas_fe.pdf', width = 16, height = 9)
embed_fonts('../Apps/ShareLaTeX/Bond reversals/betas_fe.pdf')
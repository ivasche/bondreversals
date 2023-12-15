# ==================================================================================================
# [6a] COMPUTE AND ANALYZE REVERSAL ESTMATION
#  - RUN TIME: 
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer); library(lfe); library(plm); library(extrafont)
library(sandwich); library(lmtest); library(car); library(plotrix)
library(Hmisc); library(estimatr); library(remotes)
library(ggplot2); library(cowplot); library(grid); library(ggridges); library(ggExtra)
library(scales)

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
subtrace_prefilt <- subtrace_prefilt[, .(cusip_id, ind, rating_num)]
#remove(subtrace_prefilt)
gc()

# mark corporate debentures
#cdeb <- submergent[, .(cusip_id, BOND_TYPE)][subtrace[, .N, cusip_id],, on = 'cusip_id']
#cdeb <- cdeb[, c(1,2)]
#subtrace <- cdeb[subtrace,, on = 'cusip_id']
#subtrace <- delete(subtrace, subtrace$BOND_TYPE == 'CDEB')
#remove(cdeb)
#subtrace[, BOND_TYPE := NULL]
#gc()

# print bond composition
tab = submergent[, .(cusip_id, BOND_TYPE)][subtrace[, .N, cusip_id],, on = 'cusip_id'][, .N, BOND_TYPE]
tab <- tab[order(-N)]
tab <- as.data.frame(tab[, .(100*N/sum(N))])
rownames(tab) <- c('US Corporate Debentures',
                   'US Corporate Medium-Term Notes',
                   'US Retail Notes',
                   'US Corporate Bank Note')
colnames(tab) = c('\\% of sample')

print(xtable(tab, align = c('l|', 'c'),
             digits = c(1), auto = F),
      file = '../Apps/ShareLaTeX/Bond reversals/sample_composition.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)))

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

# define function to compute stuff for equities
revfunc_eq = function(r, v) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, rv = r*v)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 12))
    
  } else {
    
    m = lm('rprime ~ r + rv', data = df)
    mod = summary(m)
    as.list(c(coef(mod)[2:3,1],
              coef(mod)[2:3,3]))
    
  }
}

# define function to compute asymmetric reversal
revfunc_asym = function(r, v, ivpos, ivneg) {
  
  ## make dataset
  df = data.frame(rprime = shift(r, type = 'lead', 1),
                  r = r, v = v, rv = r*v, rivpos = r*ivpos, rivneg = r*ivneg)
  df = df[complete.cases(df),]
  
  ## check that there's still enough data
  if ( nrow(df) < Nlim | var(df$v) == 0 ) {
    
    as.list(rep(NA_real_, 8))
    
  } else {
    
    m = lm('rprime ~ r + rv + rivpos + rivneg', data = df)
    mod = summary(m)
    as.list(c(coef(mod)[2:5,1],
              coef(mod)[2:5,3]))
    
  }
  
}


# DAILY --------------------------------------------------------------------------------------------

# apply it
subtrace[, c('C1', 'C2', 'C3', 't1', 't2', 't3', 'corv', 'pcorv',
             'delta', 'pdelta', 'corviv', 'pcorviv', 'corvivsign', 'pcorvivsign',
             'coefdiff', 'pcoefdiff') := revfunc(r = ret, v = ctrnvr_alt_3, iv = itrnvr_alt_3,
                                                 ivsign = itrnvr_alt_4),
         .(cusip_id, ind)]
gc()

# make a small summary
res = subtrace[, .(C1[1], C2[1], C3[1], t1[1], t2[1], t3[1], corv[1], pcorv[1],
                   delta[1], pdelta[1], corviv[1], pcorviv[1], corvivsign[1], pcorvivsign[1],
                   coefdiff[1], pcoefdiff[1]), .(cusip_id, issuer, ind, Nind)]
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
                    sum((C3>0)*(t3>1.65)), sum((C3<0)*(t3 < -1.65)), length(C3))])
colnames(tab) = c('~Mean~', '~Med.~', 'No.$>$0~', 'No.$<$0~', 'No.$>$0*', 'No.$<$0*', 'No. Obs.')
rownames(tab) = c('$\\hat{\\beta}_1$','$\\hat{\\beta}_2$','$\\hat{\\beta}_3$')

print(xtable(tab, align = c('l|', 'c', 'c|', 'c', 'c|', 'c', 'c|', 'c'),
             digits = c(0, 4, 4, rep(0, 5)), auto = F),
      file = '../Apps/ShareLaTeX/Bond reversals/restab.tex', floating = F,
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
      file = '../Apps/ShareLaTeX/Bond reversals/corvol.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)))

# --------------------------------------------------------------------------------------------------
# (2b) COMPUTE REVERSALS for EQUITIES
# --------------------------------------------------------------------------------------------------

# apply it
eqdf[, c('C1', 'C2', 't1', 't2') := revfunc_eq(r = RET, v = ctrnvr_alt_3), .(cusip_id, ind)]
gc()

# make a small summary
res_eq = eqdf[, lapply(.SD, mean, na.rm = T), .(cusip_id, ind), .SDcols = c(c('C1', 'C2', 't1', 't2'))]
names(res_eq)[3:6] = c('C1eq', 'C2eq', 't1eq', 't2eq')

# filter
res_eq = res_eq[!is.na(C1eq) & !is.na(C2eq)]

res_eq = res_eq[C1eq > quantile(C1eq, 0.01) & C1eq < quantile(C1eq, 0.99) &
                C2eq > quantile(C2eq, 0.01) & C2eq < quantile(C2eq, 0.99)]

# means
round(res_eq[, lapply(.SD, function(x) {coef(summary(lm(x ~ 1)))[c(1, 4)]}),,
          .SDcols = c('C1eq', 'C2eq')], 3)

# medians
round(res_eq[, lapply(.SD, median),, .SDcols = c('C1eq', 'C2eq')], 3)

# --------------------------------------------------------------------------------------------------
# (2c) COMPUTE ASYMMETRIC REVERSAL
# --------------------------------------------------------------------------------------------------

# split reversal into positive and neagtive parts
subtrace[, ivol_pos := 0]
subtrace[ivol > 0, ivol_pos := volb-vols]
subtrace[, ivol_pos := (abs(ivol_pos) - mean(abs(ivol_pos)))/sd(abs(ivol_pos)), .(cusip_id, ind)]

subtrace[, ivol_neg := 0]
subtrace[ivol < 0, ivol_neg := -(volb-vols)]
subtrace[, ivol_neg := (abs(ivol_neg) - mean(abs(ivol_neg)))/sd(abs(ivol_neg)), .(cusip_id, ind)]

subtrace[, ivol_pos := winsorize(ivol_pos, 0, 0.98)]
subtrace[, ivol_neg := winsorize(ivol_neg, 0, 0.98)]

# apply it
subtrace[, c('C1_asm', 'C2_asm', 'C3_asm', 'C4_asm',
             't1_asm', 't2_asm', 't3_asm', 't4_asm') := revfunc_asym(r = ret,
                                                                     v = ctrnvr_alt_3,
                                                                     ivpos = ivol_pos,
                                                                     ivneg = ivol_neg),
     .(cusip_id, ind)]
gc()

# make a small summary
res_asm = subtrace[, .(C1_asm[1], C2_asm[1], C3_asm[1], C4_asm[1],
                       t1_asm[1], t2_asm[1], t3_asm[1], t4_asm[1]),
                   .(cusip_id, issuer, ind, Nind)]
names(res_asm)[5:12] = c('C1_asm', 'C2_asm', 'C3_asm', 'C4_asm',
                         't1_asm', 't2_asm', 't3_asm', 't4_asm')

# filter
res_asm = res_asm[!is.na(C1_asm) & !is.na(C2_asm) & !is.na(C3_asm) & !is.na(C4_asm)]

res_asm = res_asm[C1_asm > quantile(C1_asm, 0.01) & C1_asm < quantile(C1_asm, 0.99) &
                    C2_asm > quantile(C2_asm, 0.01) & C2_asm < quantile(C2_asm, 0.99) &
                    C3_asm > quantile(C3_asm, 0.01) & C3_asm < quantile(C3_asm, 0.99) &
                    C4_asm > quantile(C4_asm, 0.01) & C4_asm < quantile(C4_asm, 0.99)]


# --------------------------------------------------------------------------------------------------
# (3) ANALYZE REVERSALS
# --------------------------------------------------------------------------------------------------

# calculate average values of bond and issuer characteristics per active period
toadd = subtrace[, lapply(.SD, function(z) {mean(z, na.rm =T)}), .(cusip_id, ind),
                 .SDcols = c('size_hist', 'rating_num', 'stock_ba', 'SE', 'VE',
                             'dealers', 'nfunds', 'avg_ba', 'avg_ba_lrg')]
res = toadd[res,,on = c('cusip_id', 'ind')]
res = res_eq[, .(cusip_id, ind, C1eq, C2eq)][res,, on = c('cusip_id', 'ind')]
res = res_asm[, .(cusip_id, ind, C1_asm, C2_asm, C3_asm, C4_asm)][res,, on = c('cusip_id', 'ind')]

# create a cross-section of bonds by averaging across different active periods for each bond
# weighting with number of trading days in each active period
csres = res[, lapply(.SD, function(x) {weighted.mean(x, w = Nind, na.rm = T)}), .(cusip_id, issuer),
            .SDcols = c('C1', 'C2', 'C3', 'corv', 'delta', 'corviv', 'C1eq', 'C2eq',
                        'C1_asm', 'C2_asm', 'C3_asm', 'C4_asm',
                        'size_hist', 'rating_num', 'stock_ba', 'SE', 'VE',
                        'dealers', 'nfunds', 'avg_ba', 'avg_ba_lrg')]

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
csres[, avg_ba_lrg_neg := -avg_ba_lrg]
csres[, stock_ba_neg := -stock_ba]
csres[, nfunds_neg := -nfunds]
csres[, dealers_neg := -dealers]
csres[, VE_neg := -VE]
csres[, size_hist_neg := -size_hist]

# PLOT scatter between bond and equity

toplot <- csres[, lapply(.SD, mean, na.rm = T), issuer, .SDcols = c('C1', 'C2', 'C3', 'C1eq', 'C2eq')]

p1 <- toplot[!is.na(C1eq)] %>%
  ggplot(aes(x = C1, y = C1eq)) + 
  geom_point() + geom_smooth(method = 'lm') +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + 
  theme_bw(base_size = 22) + 
  theme(panel.grid.major = element_line(linetype = 2, size = 0.15),
        panel.grid.minor = element_line(linetype = 2, size = 0.15),
        plot.margin = unit(c(1,14,1,1), "pt")) +
  labs(x = expression(hat(beta)[1]), y = expression(hat(beta)[1]^{eq}))

p2 <- toplot[!is.na(C1eq)] %>%
  ggplot(aes(x = C2, y = C1eq)) + 
  geom_point() + geom_smooth(method = 'lm') +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + 
  theme_bw(base_size = 22) + 
  theme(panel.grid.major = element_line(linetype = 2, size = 0.15),
        panel.grid.minor = element_line(linetype = 2, size = 0.15),
        plot.margin = unit(c(1,14,1,1), "pt")) +
  labs(x = expression(hat(beta)[2]), y = expression(hat(beta)[1]^{eq}))

p3 <- toplot[!is.na(C1eq)] %>%
  ggplot(aes(x = C3, y = C1eq)) + 
  geom_point() + geom_smooth(method = 'lm') +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + 
  theme_bw(base_size = 22) + 
  theme(panel.grid.major = element_line(linetype = 2, size = 0.15),
        panel.grid.minor = element_line(linetype = 2, size = 0.15),
        plot.margin = unit(c(1,14,1,1), "pt")) +
  labs(x = expression(hat(beta)[3]), y = expression(hat(beta)[1]^{eq}))

p4 <- toplot[!is.na(C1eq)] %>%
  ggplot(aes(x = C1, y = C2eq)) + 
  geom_point() + geom_smooth(method = 'lm') +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + 
  theme_bw(base_size = 22) + 
  theme(panel.grid.major = element_line(linetype = 2, size = 0.15),
        panel.grid.minor = element_line(linetype = 2, size = 0.15),
        plot.margin = unit(c(1,14,1,1), "pt")) +
  labs(x = expression(hat(beta)[1]), y = expression(hat(beta)[2]^{eq}))

p5 <- toplot[!is.na(C1eq)] %>%
  ggplot(aes(x = C2, y = C2eq)) + 
  geom_point() + geom_smooth(method = 'lm') +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + 
  theme_bw(base_size = 22) + 
  theme(panel.grid.major = element_line(linetype = 2, size = 0.15),
        panel.grid.minor = element_line(linetype = 2, size = 0.15),
        plot.margin = unit(c(1,14,1,1), "pt")) +
  labs(x = expression(hat(beta)[2]), y = expression(hat(beta)[2]^{eq}))

p6 <- toplot[!is.na(C1eq)] %>%
  ggplot(aes(x = C3, y = C2eq)) + 
  geom_point() + geom_smooth(method = 'lm') +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + 
  theme_bw(base_size = 22) + 
  theme(panel.grid.major = element_line(linetype = 2, size = 0.15),
        panel.grid.minor = element_line(linetype = 2, size = 0.15),
        plot.margin = unit(c(1,14,1,1), "pt")) +
  labs(x = expression(hat(beta)[3]), y = expression(hat(beta)[2]^{eq}))

namef = '../Apps/ShareLaTeX/Bond reversals/betahats_bond_equity.pdf'
pdf(file = namef, family = 'CM Roman', width = 16, height = 9, pointsize = 24)
  plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)
dev.off()
embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))

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
pc_cols_4 = c('VE_neg', 'stock_ba')

pc_cols = list(pc_cols_1, pc_cols_2, pc_cols_3, pc_cols_4)
pcs = lapply(pc_cols, function(x) {pc_func(csdata = csres, vrbls = x)})
names(pcs) = paste0('pc', 1:length(pcs))
pcs = as.data.table(pcs)

csres = cbind(csres, pcs)

# some investigation of PCs

pc_test = lapply(pc_cols[1:3], function(x) {
    prcomp(formula = formula(paste('~', paste(x, collapse = '+'))),
           data = csres, center = T, scale = T, na.action = na.exclude)
  }) 

tab = t(rbind(pc_test[[1]]$rotation[, 1],
            c(pc_test[[2]]$rotation[, 1], rep(NA_real_, 2)),
            c(NA_real_, (-1)*pc_test[[3]]$rotation[, 1], rep(NA_real_, 2))
))
tab = rbind(tab, sapply(pc_test, function(x) {x$sdev[1]^2/sum(x$sdev^2)}))  
rownames(tab) = c('Bond bid-ask',
                  '--No. mutual fund owners ',
                  '--No. dealers',
                  '--Issue size',
                  '--Issuer size',
                  'Stock bid-ask',
                  'Share of explained variance, \\%'
                  )
colnames(tab) = c('~~~PC$_\\text{all}~~~~$',
                  '~~~PC$_\\text{bond}~~~$',
                  'PC$_\\text{bond-ex-ba}$')

print(xtable(tab, align = c('l|', 'c', 'c', 'c'),
             digits = 2),
      file = '../Apps/ShareLaTeX/Bond reversals/pcs.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)-1, nrow(tab)))

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
             auto = F), file = '../Apps/ShareLaTeX/Bond reversals/sumstat_cs_D.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, 3, 4, 10, nrow(tab)))

tab = corstars(csres[, .(avg_ba, nfunds, size_hist, dealers, VE, stock_ba, pc1, pc2, pc3)])
colnames(tab) =  c('Bond bid-ask', 'No. funds', 'Issue size', 'No. dealers', 'Issuer size',
                   'Stock bid-ask', 'PC\\textsubscript{all}', 'PC\\textsubscript{bond}')
rownames(tab) =  c('Bond bid-ask', 'No. funds', 'Issue size', 'No. dealers', 'Issuer size', 'Stock bid-ask',
                   'PC\\textsubscript{all}', 'PC\\textsubscript{bond}', 'PC\\textsubscript{bond-ex-ba}')

print(xtable(tab[2:9,], align = c('l|', rep('c', 8)),
             auto = F), file = '../Apps/ShareLaTeX/Bond reversals/sumcorr_cs_D.tex', floating = F,
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
  
  fmla[[1]] = formula(paste(coef, '~ avg_ba | rating | 0 | rating'))
  fmla[[2]] = formula(paste(coef, '~ nfunds_neg | rating | 0 | rating'))
  fmla[[3]] = formula(paste(coef, '~ size_hist_neg | rating | 0 | rating'))
  fmla[[4]] = formula(paste(coef, '~ dealers_neg | rating | 0 | rating'))
  fmla[[5]] = formula(paste(coef, '~ VE_neg | rating | 0 | rating'))
  fmla[[6]] = formula(paste(coef, '~ stock_ba | rating | 0 | rating'))

  fmla[[7]] = formula(paste(coef, '~ avg_ba + nfunds_neg + size_hist_neg + dealers_neg | rating | 0 | rating'))
  fmla[[8]] = formula(paste(coef, '~ avg_ba + nfunds_neg + size_hist_neg + dealers_neg + VE_neg + stock_ba | rating | 0 | rating'))
  
  fmla[[9]] = formula(paste(coef, '~ pc1 | rating | 0 | rating'))
  fmla[[10]] = formula(paste(coef, '~ pc2 | rating | 0 | rating'))
  fmla[[11]] = formula(paste(coef, '~ pc3 | rating | 0 | rating'))

  fmla[[12]] = formula(paste(coef, '~ corv + corviv + pc1 | rating | 0 | rating'))
  fmla[[13]] = formula(paste(coef, '~ corv + corviv + pc2 | rating | 0 | rating'))
  fmla[[14]] = formula(paste(coef, '~ corv + corviv + pc3 | rating | 0 | rating'))
  
  fmla[[15]] = formula(paste(coef, '~ pc4 | rating | 0 | rating'))

  mlist = lapply(fmla, function(x) {felm(x, data = as.data.frame(csres), keepX = T,
                                         cmethod='reghdfe')})
}

# run estimations 
est_c1 = estfunc(coef = 'C1')
est_c2 = estfunc(coef = 'C2')
est_c3 = estfunc(coef = 'C3')
est_c1eq = estfunc(coef = 'C1eq')
est_c2eq = estfunc(coef = 'C2eq')
est_c1asm = estfunc(coef = 'C1_asm')
est_c2asm = estfunc(coef = 'C2_asm')
est_c3asm = estfunc(coef = 'C3_asm')
est_c4asm = estfunc(coef = 'C4_asm')

# save output for c1
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test.tex', est_c1[1:8],
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
            add.lines = list(c("Rating FE", rep('YES', length(est_c1[1:8])))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test.tex')

# save output for c2
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test2.tex', est_c2[1:8],
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
            add.lines = list(c("Rating FE", rep('YES', length(est_c1[1:8])))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test2.tex')
  
  
# save output for c3
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test3.tex', est_c3[1:8],
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
            add.lines = list(c("Rating FE", rep('YES', length(est_c1[1:8])))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test3.tex')

# save joint
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc.tex',
          c(est_c1[9:11], est_c2[9:11], est_c3[9:11]),
          font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
          covariate.labels = c('PC$_\\text{all}$',
                               'PC$_\\text{bond}$',
                               'PC$_\\text{bond-ex-ba}$'),
          float = F, intercept.top = T, intercept.bottom = F, digits = 3,
          dep.var.labels.include = F, model.numbers = T,
          column.labels = c('$\\hat{\\beta}_1$' , '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$'),
          column.separate = rep(3,3),          
          dep.var.caption = '',
          add.lines = list(c("Rating FE", rep('YES', ncol(pcs)*3))),
          notes.align = 'l')
)
mod <- sub("lccccccccc", "lccc|ccc|ccc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc.tex')

# save output for equities
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_equity.tex',
            c(est_c1eq[c(5,6,15)], est_c2eq[c(5,6,15)]),
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('$\\text{--Issuer size}$',
                                 '$\\text{Stock bid-ask}$',
                                 'PC$_\\text{equity}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = F, model.numbers = T,
            column.labels = c('$\\hat{\\beta}_1^{eq}$' , '$\\hat{\\beta}_2^{eq}$'),
            column.separate = c(3,3),          
            dep.var.caption = '',
            add.lines = list(c("Rating FE", rep('YES', 6))),
            notes.align = 'l')
)
mod <- sub("lcccccc", "lccc|ccc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_equity.tex')

# save output for asymmetric beta3-beta4
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_asym.tex',
            c(est_c1asm[c(10)], est_c2asm[c(10)], est_c3asm[c(10)], est_c4asm[c(10)]),
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = T, model.numbers = F,
            dep.var.labels = c('$\\hat{\\beta}$\\textsubscript 1',
                               '$\\hat{\\beta}$\\textsubscript 2',
                               '$\\hat{\\beta}$\\textsubscript {3,+}',
                               '$\\hat{\\beta}$\\textsubscript {3,--}'),
            dep.var.caption = '',
            add.lines = list(c("Rating FE", rep('YES', 4))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_asym.tex')


# save joint with extra controls
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc_vlm.tex',
          c(est_c1[12:14], est_c2[12:14], est_c3[12:14]),
          font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
          covariate.labels = c('$\\text{C-to-C vlm. corr.}$',
                               '$\\text{C-to-D vlm. corr.}$',
                               'PC$_\\text{all}$',
                               'PC$_\\text{bond}$',
                               'PC$_\\text{bond-ex-ba}$'),
          float = F, intercept.top = T, intercept.bottom = F, digits = 3,
          dep.var.labels.include = F, model.numbers = F,
          column.labels = c('$\\hat{\\beta}_1$' , '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$'),
          column.separate = rep(3,3),          
          dep.var.caption = '',
          add.lines = list(c("Rating FE", rep('YES', ncol(pcs)*3))),
          notes.align = 'l')
)
mod <- sub("lccccccccc", "lccc|ccc|ccc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc_vlm.tex')

vlm_robust <- sapply(as.list(c(est_c1[13], est_c2[13], est_c3[13],
                               est_c1[14], est_c2[14], est_c3[14])),
                     function(x) {summary(x)$coefficients[3, c(1,2,4)]}
)
save(vlm_robust, file = './data_input/vlm_robust.RData')

# same for industries
toadd <- submergent[, .(cusip_id, INDUSTRY_CODE)]
csres <- toadd[csres,, on = 'cusip_id']
csres[INDUSTRY_CODE %in% c(seq(10, 16, 1), 32), industry := 'Industrial']
csres[INDUSTRY_CODE %in% seq(20, 26, 1), industry := 'Financial']
csres[INDUSTRY_CODE %in% c(30, 31, 33), industry := 'Utility']

estfunc_industry = function(coef) {
  
  fmla = formula(paste(coef, '~ pc2 | rating | 0 | rating'))

  mlist = lapply(as.list(c('Industrial', 'Financial', 'Utility')),
                 function(x) {felm(formula = fmla,
                                   data = as.data.frame(csres[industry == x]),
                                   keepX = T,
                                   cmethod='reghdfe')})
  names(mlist) <- c('Industrial', 'Financial', 'Utility')
  return(mlist)
  
}

est_c1_industry = estfunc_industry(coef = 'C1')
est_c2_industry = estfunc_industry(coef = 'C2')
est_c3_industry = estfunc_industry(coef = 'C3')

# save industry
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_industry.tex',
            c(est_c1_industry[1], est_c2_industry[1], est_c3_industry[1],
              est_c1_industry[2], est_c2_industry[2], est_c3_industry[2],
              est_c1_industry[3], est_c2_industry[3], est_c3_industry[3]),
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = T, model.numbers = F,
            dep.var.labels = rep(c('$\\hat{\\beta}$\\textsubscript 1',
                                   '$\\hat{\\beta}$\\textsubscript 2',
                                   '$\\hat{\\beta}$\\textsubscript 3'),
                                 3),
            column.labels = c('Industrial', 'Financial', 'Utility'),
            column.separate = rep(3,3,3),          
            dep.var.caption = '',
            add.lines = list(c("Rating FE", rep('YES', 9))),
            notes.align = 'l')
)
mod <- sub("lccccccccc", "lccc|ccc|ccc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_industry.tex')


# run the same model with lm_robust and rating fixed effects (for plotting)

estfunc_lmrob = function(coef) {
  
  fmla = vector('list', 3)

  fmla[[1]] = formula(paste(coef, '~ pc1 + factor(rating) - 1'))
  fmla[[2]] = formula(paste(coef, '~ pc2 + factor(rating) - 1'))
  fmla[[3]] = formula(paste(coef, '~ pc3 + factor(rating) - 1'))
  
  mlist = lapply(fmla, function(x) {lm_robust(formula = x, data = csres, clusters = rating)})
}

est_c1_lmrob = estfunc_lmrob(coef = 'C1')
est_c2_lmrob = estfunc_lmrob(coef = 'C2')
est_c3_lmrob = estfunc_lmrob(coef = 'C3')

# predict from the model with pc2

pcqs = quantile(csres[, pc2], seq(0.10, 0.90, 0.05), na.rm = T)
#pcqs = pcqs[(length(pcqs):1)]

revcomp = vector('list', 5)

# b1
pb1 = predict(object = est_c1_lmrob[[2]], interval = 'confidence', data.frame(pc2 = pcqs, rating = 8))
revcomp$b1 <- t(as.matrix(pb1$fit))

# b2
pb2 = predict(object = est_c2_lmrob[[2]], interval = 'confidence', data.frame(pc2 = pcqs, rating = 8))
revcomp$b2 <- t(as.matrix(pb2$fit))

# b3
pb3 = predict(object = est_c3_lmrob[[2]], interval = 'confidence', data.frame(pc2 = pcqs, rating = 8))
revcomp$b3 <- t(as.matrix(pb3$fit))

volmult = 2

revcomp$highcvol  = revcomp$b1 + volmult*revcomp$b2
revcomp$highivol  = revcomp$b1 + volmult*revcomp$b3

# plot stuff with betas
range = 1:length(pcqs)

namef = '../Apps/ShareLaTeX/Bond reversals/betahat.pdf'
pdf(file = namef, family = 'CM Roman', width = 16, height = 6, pointsize = 26)

par(mfrow = c(1, 3), mar = c(5,4,2.5,0.5))

plot(x = range, y = revcomp$b1[1,range], type = 'n',
     ylim = c(min(revcomp$b1[, range]), max(revcomp$b1[, range])),
     xlab = '', ylab = '', las = 1, xaxt = 'n',
     main = bquote(paste('E[', hat(beta[1]), '|covariates]')))
abline(h = seq(-0.40, -0.25, 0.05), lty = 3, col = rgb(0,0,0,0.25))
lines(x = range, y = revcomp$b1[1,range], lwd = 3, col = 'black')
polygon(x = c(range,rev(range)), y = c(revcomp$b1[2,range], revcomp$b1[3,rev(range)]),
        col = rgb(0,0,1, 0.15), border = NA)
mtext(text = 'High', side = 1, line = 1, at = 16, cex = 0.8)
mtext(text = 'Average', side = 1, line = 1, at = 9, cex = 0.8)
mtext(text = 'Low', side = 1, line = 1, at = 2, cex = 0.8)
#mtext(text = 'Info asymmetry', side = 1, line = 3, at = 5, cex = 0.8)

plot(x = range, y = revcomp$b2[1,range], type = 'n',
     #     ylim = c(min(revcomp$b3[, range]), max(revcomp$b3[, range])),
     ylim = c(0.02, 0.10),
     xlab = '', ylab = '', las = 1, xaxt = 'n',
     main = bquote(paste('E[', hat(beta[2]), '|covariates]')))
abline(h = seq(0.02, 0.10, 0.02), lty = 3, col = rgb(0,0,0,0.25))
lines(x = range, y = revcomp$b2[1,range], lwd = 3, col = 'black')
polygon(x = c(range,rev(range)), y = c(revcomp$b2[2,range], revcomp$b2[3,rev(range)]),
        col = rgb(0,0,1, 0.15), border = NA)
mtext(text = 'High', side = 1, line = 1, at = 16, cex = 0.8)
mtext(text = 'Average', side = 1, line = 1, at = 9, cex = 0.8)
mtext(text = 'Low', side = 1, line = 1, at = 2, cex = 0.8)
mtext(text = 'Info asymmetry', side = 1, line = 3, at = 7, cex = 0.8)

plot(x = range, y = revcomp$b3[1,range], type = 'n',
#     ylim = c(min(revcomp$b3[, range]), max(revcomp$b3[, range])),
      ylim = c(0.02, 0.10),
      xlab = '', ylab = '', las = 1, xaxt = 'n',
      main = bquote(paste('E[', hat(beta[3]), '|covariates]')))
abline(h = seq(0.02, 0.10, 0.02), lty = 3, col = rgb(0,0,0,0.25))
lines(x = range, y = revcomp$b3[1,range], lwd = 3, col = 'black')
polygon(x = c(range,rev(range)), y = c(revcomp$b3[2,range], revcomp$b3[3,rev(range)]),
        col = rgb(0,0,1, 0.15), border = NA)
mtext(text = 'High', side = 1, line = 1, at = 16, cex = 0.8)
mtext(text = 'Average', side = 1, line = 1, at = 9, cex = 0.8)
mtext(text = 'Low', side = 1, line = 1, at = 2, cex = 0.8)
#mtext(text = 'Info asymmetry', side = 1, line = 3, at = 5, cex = 0.8)

dev.off()
embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))

# plot stuff with volumes

namef = '../Apps/ShareLaTeX/Bond reversals/volumes.pdf'
pdf(file = namef, family = 'CM Roman', width = 16, height = 6, pointsize = 26)

par(mfrow = c(1, 3), mar = c(5,4,2.5,0.5))

plot(x = range, y = revcomp$b1[1,range], type = 'n',
     ylim = c(-0.45, -0.10),
     xlab = '', ylab = '', las = 1, xaxt = 'n',
     main = 'Average volume day')
abline(h = seq(-0.45, -0.15, 0.05), lty = 3, col = rgb(0,0,0,0.25))
lines(x = range, y = revcomp$b1[1,range], lwd = 3, col = 'black')
#polygon(x = c(range,rev(range)), y = c(revcomp$b1[2,range], revcomp$b1[3,rev(range)]),
#        col = rgb(0,0,1, 0.15), border = NA)
mtext(text = 'High', side = 1, line = 1, at = 16, cex = 0.8)
mtext(text = 'Average', side = 1, line = 1, at = 9, cex = 0.8)
mtext(text = 'Low', side = 1, line = 1, at = 2, cex = 0.8)
#mtext(text = 'Info asymmetry', side = 1, line = 3, at = 5, cex = 0.8)

plot(x = range, y = revcomp$highcvol[1,range], type = 'n',
     ylim = c(-0.45, -0.10), 
     xlab = '', ylab = '', las = 1, xaxt = 'n',
     main = 'High C-to-C volume day')
abline(h = seq(-0.45, -0.15, 0.05), lty = 3, col = rgb(0,0,0,0.25))
lines(x = range, y = revcomp$highcvol[1,range], lwd = 3, col = 'black')
#polygon(x = c(range,rev(range)),
#        y = c(revcomp$highcvol[2,range], revcomp$highcvol[3,rev(range)]),
#        col = rgb(0,0,1, 0.15), border = NA)
mtext(text = 'High', side = 1, line = 1, at = 16, cex = 0.8)
mtext(text = 'Average', side = 1, line = 1, at = 9, cex = 0.8)
mtext(text = 'Low', side = 1, line = 1, at = 2, cex = 0.8)
mtext(text = 'Info asymmetry', side = 1, line = 3, at = 7, cex = 0.8)

plot(x = range, y = revcomp$highivol[1,range], type = 'n',
     ylim = c(-0.45, -0.10), 
     xlab = '', ylab = '', las = 1, xaxt = 'n',
     main = 'High C-to-D volume day')
abline(h = seq(-0.45, -0.15, 0.05), lty = 3, col = rgb(0,0,0,0.25))
lines(x = range, y = revcomp$highivol[1,range], lwd = 3, col = 'black')
#polygon(x = c(range,rev(range)),
#        y = c(revcomp$highivol[2,range], revcomp$highivol[3,rev(range)]),
#        col = rgb(0,0,1, 0.15), border = NA)
mtext(text = 'High', side = 1, line = 1, at = 16, cex = 0.8)
mtext(text = 'Average', side = 1, line = 1, at = 9, cex = 0.8)
mtext(text = 'Low', side = 1, line = 1, at = 2, cex = 0.8)
#mtext(text = 'Info asymmetry', side = 1, line = 3, at = 5, cex = 0.8)

dev.off()
embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))

# evidence for big issuers -------------------------------------------------------------------------

estfunc_big = function(coef) {
  
  df = csres
  df = df[issuer %in% df[,.(.N), issuer][N>15, issuer]]
  
  fmla = vector('list', 3)
  
  fmla[[1]] = formula(paste(coef, '~ pc2 | rating | 0 | rating'))
  fmla[[2]] = formula(paste(coef, '~ pc2 | issuer + rating | 0 | rating'))
  fmla[[3]] = formula(paste(coef, '~ pc2 | issuer + rating | 0 | issuer'))

  mlist = lapply(fmla, function(x) {felm(x, data = as.data.frame(df), keepX = T,
                                         cmethod='reghdfe')})
}

# run estimations 
est_c1_big = estfunc_big(coef = 'C1')
est_c2_big = estfunc_big(coef = 'C2')
est_c3_big = estfunc_big(coef = 'C3')

# save output
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_pc_big.tex',
            c(est_c1_big[3], est_c2_big[3], est_c3_big[3]),
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('\\hline PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = F, model.numbers = F,
            column.labels = c('$\\hat{\\beta}$\\textsubscript 1',
                              '$\\hat{\\beta}$\\textsubscript 2',
                              '$\\hat{\\beta}$\\textsubscript 3'),
            dep.var.caption = '',
            add.lines = list(c("Rating, Issuer FE", rep('YES', 3)),
                             c("Issuer-clustered SE", rep(c('YES'), 3))),
            notes.align = 'l')
)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_pc_big.tex')


# evidence for HY and IG separately ----------------------------------------------------------------

estfunc_ighy = function(coef) {
  
  df_ig = csres[rating <= 10]
  df_hy = csres[rating > 10]
  
  fmla = vector('list')
  
  fmla[[1]] = formula(paste(coef, '~ pc2 | rating | 0 | rating'))

  mlist = c(lapply(fmla, function(x) {felm(x, data = as.data.frame(df_ig), keepX = T,
                                         cmethod='reghdfe')}),
            lapply(fmla, function(x) {felm(x, data = as.data.frame(df_hy), keepX = T,
                                           cmethod='reghdfe')}))

}

# run estimations 
est_c1_ighy = estfunc_ighy(coef = 'C1')
est_c2_ighy = estfunc_ighy(coef = 'C2')
est_c3_ighy = estfunc_ighy(coef = 'C3')

# save output
mod <- capture.output(
  stargazer(out = '../Apps/ShareLaTeX/Bond reversals/test_ighy.tex',
            c(est_c1_ighy, est_c2_ighy, est_c3_ighy),
            font.size = "scriptsize", no.space = T, omit.stat = c('ser', 'f', 'adj.rsq'),
            covariate.labels = c('PC$_\\text{bond}$'),
            float = F, intercept.top = T, intercept.bottom = F, digits = 3,
            dep.var.labels.include = T, model.numbers = F,
            column.labels = rep(c('IG', 'HY'), 3),
            dep.var.caption = '',
            dep.var.labels = c('$\\hat{\\beta}_1$', '$\\hat{\\beta}_2$', '$\\hat{\\beta}_3$'),
            add.lines = list(c("Rating FE", rep(rep('YES', length(est_c1_big)), 3))),
            notes.align = 'l')
)
mod <- sub("lcccccc", "lcc|cc|cc", mod)
mod <- sub("5pt", "0pt", mod)
writeLines(mod, '../Apps/ShareLaTeX/Bond reversals/test_ighy.tex')

# --------------------------------------------------------------------------------------------------
# (4) STRATEGY
# --------------------------------------------------------------------------------------------------

# add transaction costs
load('./data_input/liquidity_irtc.RData')
irtc[, trd_exctn_dt := as.IDate(as.Date(as.character(trd_exctn_dt), format = '%Y%m%d'))]

irtc <- irtc[, weighted.mean(x = irtc/100, w = nobs), .(cusip_id, as.yearmon(trd_exctn_dt))]
names(irtc)[3] <- 'irtc'
irtc[irtc == 0, irtc := NA]
names(irtc)[2] <- 'ym'
gc()

# take a monthly dataset for performance measurement
submex_init[, retm := winsorize(retm, 0.001, 0.999)]
submex = submex_init[ym >= as.yearmon('2005-01') & ym <= as.yearmon('2018-12')]

# add
submex <- irtc[submex,, on = c('cusip_id', 'ym')]
remove(irtc)
gc()

# make index variable
#submex[is.na(index) | index < 0.5, index := 0]
#submex[is.na(index) | index >= 0.5, index := 1]

# make lagged variables
setkey(submex, cusip_id, ym)
submex[, lretm := shift(retm), cusip_id]
submex[, lsize := shift(size), cusip_id]
submex[, lrat  := shift(rat), cusip_id]
#submex[, lindex  := shift(index, type = 'lag', n = 4L), cusip_id]
submex[, lnfunds  := shift(nfunds, type = 'lag', n = 7L), cusip_id]
submex[, lcvol     := shift(cvol), cusip_id]
submex[, lVE     := shift(VE), cusip_id]

# 12-month backward looking avg ba
submex[, ba := (if (.N < 12) {
  rep(NA_real_, .N)
} else {
  rollapply(avg_ba/100, function(x) {mean(x, na.rm = T)}, width = 12, fill = NA, align = 'right')
}), cusip_id]
submex[, lba := shift(ba), cusip_id]

# 12-month backward looking roundtrip
submex[, ba_alt := (if (.N < 12) {
  rep(NA_real_, .N)
} else {
  rollapply(irtc/100, function(x) {mean(x, na.rm = T)}, width = 12, fill = NA, align = 'right')
}), cusip_id]
submex[, lba_alt := shift(ba_alt), cusip_id]

# save before filtering
submex_prefilt = submex
submex = submex_prefilt

# remove na lags
submex = submex[!is.na(lsize) & !is.na(lnfunds) & lba < 0.01 & lrat <= 21 & lsize >= 200]

# load a risk-free rate
rf = fread('./data_input/GS3M.csv')
rf[, DATE := as.yearmon(DATE)]
rf[, GS3M := GS3M/1200]
names(rf) = c('ym', 'rf')
rf = as.xts(rf[, rf], order.by = rf[, ym])

# define quantiles
submex[, IG := 0]
submex[lrat <= 10, IG := 1]

submex[, Qrev    := .bincode(-lretm, breaks = quantile(-lretm,seq(0, 1, 0.2), na.rm = T),
                                      include.lowest=T), ym]
submex[, Qrat    := .bincode(lrat, breaks = quantile(lrat,seq(0, 1, 1/3), na.rm = T),
                                      include.lowest=T), ym]

submex[IG == 1, Qrat_IG    := .bincode(lrat, breaks = quantile(lrat,seq(0, 1, 1/3), na.rm = T),
                             include.lowest=T), ym]
submex[IG == 0, Qrat_HY    := .bincode(lrat, breaks = quantile(lrat,seq(0, 1, 1/3), na.rm = T),
                                include.lowest=T), ym]

submex[, Qnfunds := .bincode(lnfunds, breaks = quantile(lnfunds,seq(0, 1, 0.5), na.rm = T),
                                      include.lowest=T), ym]
submex[, Qnfunds := Qnfunds - 1]

submex[, QVE := .bincode(lVE, breaks = quantile(lVE,seq(0, 1, 0.5), na.rm = T),
                             include.lowest=T), ym]
submex[, QVE := QVE - 1]

# Compute returns of sorted rating and reversal portfolios conditioning on number of fund owners
# s.rat.rev = submex[!is.na(Qrev) & !is.na(Qrat),
#                   .(weighted.mean(x = retm, w = lsize, na.rm = T),
#                     weighted.mean(x = retm*Qnfunds, w = lsize*Qnfunds, na.rm = T),
#                     weighted.mean(x = retm*(1-Qnfunds), w = lsize*(1-Qnfunds), na.rm = T)),
#                   by = .(ym, Qrat, Qrev)][order(ym)]

# names(s.rat.rev)[4:6] = c('total', 'manyfunds', 'fewfunds')

# REV long
#REV = s.rat.rev[Qrev==5, lapply(.SD, mean, na.rm = T), .(ym),
#                .SDcols = c('total', 'manyfunds', 'fewfunds')]

# REV xts
#REVxts = as.xts(as.data.frame(REV[, -1]), order.by = REV[, ym])
#REVxts = merge(REVxts, mkt[, mkt])
#names(REVxts)[ncol(REVxts)] = 'mkt'
#REVxtscum = do.call(merge, lapply(REVxts+1, cumprod))

# portfolio composition
comp = vector('list', 15)
names(comp) = c(paste0('R1D_', c('all', 'many', 'few', 'big', 'small', 'IG', 'HY')),
                paste0('R2D_', c('all', 'many', 'few', 'big', 'small', 'IG', 'HY')),
                'mkt')

## R1D_all
comp$R1D_all = split(rbindlist(list(submex[Qrev == 5,
                                             .(cusip_id, lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R1D_all = comp$R1D_all[order(as.yearmon(names(comp$R1D_all)))]
for (i in 1:length(comp$R1D_all)) {names(comp$R1D_all[[i]])[3] = 'weight'}

## R1D_IG
comp$R1D_IG = split(rbindlist(list(submex[Qrev == 5 & IG == 1,
                                           .(cusip_id, lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R1D_IG = comp$R1D_IG[order(as.yearmon(names(comp$R1D_IG)))]
for (i in 1:length(comp$R1D_IG)) {names(comp$R1D_IG[[i]])[3] = 'weight'}

## R1D_HY
comp$R1D_HY = split(rbindlist(list(submex[Qrev == 5 & IG == 0,
                                          .(cusip_id, lsize/sum(lsize), retm), ym])),
                    by = 'ym')
comp$R1D_HY = comp$R1D_HY[order(as.yearmon(names(comp$R1D_HY)))]
for (i in 1:length(comp$R1D_HY)) {names(comp$R1D_HY[[i]])[3] = 'weight'}

## R1D_many
comp$R1D_many = split(rbindlist(list(submex[Qrev == 5 & Qnfunds == 1,
                                          .(cusip_id, lsize/sum(lsize), retm), ym])),
                    by = 'ym')
comp$R1D_many = comp$R1D_many[order(as.yearmon(names(comp$R1D_many)))]
for (i in 1:length(comp$R1D_many)) {names(comp$R1D_many[[i]])[3] = 'weight'}

## R1D_few
comp$R1D_few = split(rbindlist(list(submex[Qrev == 5 & Qnfunds == 0,
                                            .(cusip_id, lsize/sum(lsize), retm), ym])),
                      by = 'ym')
comp$R1D_few = comp$R1D_few[order(as.yearmon(names(comp$R1D_few)))]
for (i in 1:length(comp$R1D_few)) {names(comp$R1D_few[[i]])[3] = 'weight'}

## R1D_big
comp$R1D_big = split(rbindlist(list(submex[Qrev == 5 & QVE == 1,
                                           .(cusip_id, lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R1D_big = comp$R1D_big[order(as.yearmon(names(comp$R1D_big)))]
for (i in 1:length(comp$R1D_big)) {names(comp$R1D_big[[i]])[3] = 'weight'}
## R1D_small
comp$R1D_small = split(rbindlist(list(submex[Qrev == 5 & QVE == 0,
                                           .(cusip_id, lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R1D_small = comp$R1D_small[order(as.yearmon(names(comp$R1D_small)))]
for (i in 1:length(comp$R1D_small)) {names(comp$R1D_small[[i]])[3] = 'weight'}

## R2D_all
comp$R2D_all = split(rbindlist(list(submex[Qrev == 5 & Qrat == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat == 2,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat == 3,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R2D_all = comp$R2D_all[order(as.yearmon(names(comp$R2D_all)))]
for (i in 1:length(comp$R2D_all)) {names(comp$R2D_all[[i]])[3] = 'weight'}

## R2D_IG
comp$R2D_IG = split(rbindlist(list(submex[Qrev == 5 & Qrat_IG == 1 & IG == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat_IG == 2 & IG == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat_IG == 3 & IG == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R2D_IG = comp$R2D_IG[order(as.yearmon(names(comp$R2D_IG)))]
for (i in 1:length(comp$R2D_IG)) {names(comp$R2D_IG[[i]])[3] = 'weight'}

## R2D_HY
comp$R2D_HY = split(rbindlist(list(submex[Qrev == 5 & Qrat_HY == 1 & IG == 0,
                                          .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                   submex[Qrev == 5 & Qrat_HY == 2 & IG == 0,
                                          .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                   submex[Qrev == 5 & Qrat_HY == 3 & IG == 0,
                                          .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                    by = 'ym')
comp$R2D_HY = comp$R2D_HY[order(as.yearmon(names(comp$R2D_HY)))]
for (i in 1:length(comp$R2D_HY)) {names(comp$R2D_HY[[i]])[3] = 'weight'}

## R2D_many
comp$R2D_many = split(rbindlist(list(submex[Qrev == 5 & Qrat == 1 & Qnfunds == 1,
                                               .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                        submex[Qrev == 5 & Qrat == 2 & Qnfunds == 1,
                                               .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                        submex[Qrev == 5 & Qrat == 3 & Qnfunds == 1,
                                               .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                         by = 'ym')
comp$R2D_many = comp$R2D_many[order(as.yearmon(names(comp$R2D_many)))]
for (i in 1:length(comp$R2D_many)) {names(comp$R2D_many[[i]])[3] = 'weight'}

## R2D_few
comp$R2D_few = split(rbindlist(list(submex[Qrev == 5 & Qrat == 1 & Qnfunds == 0,
                                              .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                       submex[Qrev == 5 & Qrat == 2 & Qnfunds == 0,
                                              .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                       submex[Qrev == 5 & Qrat == 3 & Qnfunds == 0,
                                              .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                        by = 'ym')
comp$R2D_few = comp$R2D_few[order(as.yearmon(names(comp$R2D_few)))]
for (i in 1:length(comp$R2D_few)) {names(comp$R2D_few[[i]])[3] = 'weight'}

## R2D_big
comp$R2D_big = split(rbindlist(list(submex[Qrev == 5 & Qrat == 1 & QVE == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat == 2 & QVE == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat == 3 & QVE == 1,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R2D_big = comp$R2D_big[order(as.yearmon(names(comp$R2D_big)))]
for (i in 1:length(comp$R2D_big)) {names(comp$R2D_big[[i]])[3] = 'weight'}

## R2D_small
comp$R2D_small = split(rbindlist(list(submex[Qrev == 5 & Qrat == 1 & QVE == 0,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat == 2 & QVE == 0,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym],
                                    submex[Qrev == 5 & Qrat == 3 & QVE == 0,
                                           .(cusip_id, 1/3*lsize/sum(lsize), retm), ym])),
                     by = 'ym')
comp$R2D_small = comp$R2D_small[order(as.yearmon(names(comp$R2D_small)))]
for (i in 1:length(comp$R2D_small)) {names(comp$R2D_small[[i]])[3] = 'weight'}

## market
comp$mkt = split(submex[, .(cusip_id, lsize/sum(lsize), retm), ym], by = 'ym')
comp$mkt = comp$mkt[order(as.yearmon(names(comp$mkt)))]
for (i in 1:length(comp$mkt)) {names(comp$mkt[[i]])[3] = 'weight'}

# turnover and net return measurement function
# and the end of each month you rebalance: abs values of weight changes are multiplied by 
# half-spread
tcost = function(x) {
  
  seq = seq_along(names(x))
  
  A = lapply(as.list(seq), function(ind) {
    
    # assume that in the beginning of the first month you need to buy the enitre portfolio
    if (ind == seq[[1]]) {
      
      # compute rebalancing
      old = x[[ind]]; new = x[[ind]]
      names(old)[3] = 'w0'; names(new)[3] = 'w1'
      df = merge(old[,-c(1,4)], new[,-c(1,4)], all = T)
      
      # add bid-ask spreads and compute trading costs
      df = submex_prefilt[ym == as.yearmon(names(x)[ind]),
                          .(cusip_id, ba, lba, ba_alt, lba_alt, retm)][df,,on = 'cusip_id']
      df[is.na(ba), ba := median(df[, ba], na.rm = T)]
      df[is.na(lba), lba := median(df[, lba], na.rm = T)]
      df[is.na(ba_alt), ba_alt := median(df[, ba_alt], na.rm = T)]
      df[is.na(lba_alt), lba_alt := median(df[, lba_alt], na.rm = T)]
      df[, .(weighted.mean(retm, w = w1, na.rm = T),
             sum(w0),
             sum(w0*lba/2),
             sum(w0*lba_alt/2),
             nrow(df))]  
      
    } else {
      
      # compute rebalancing
      old = x[[ind-1]]; new = x[[ind]]
      names(old)[3] = 'w0'; names(new)[3] = 'w1'
      df = merge(old[, -c(1, 4)], new[, -c(1,4)], all = T)
      df[is.na(w0), w0 := 0]
      df[is.na(w1), w1 := 0]
      
      # add bid-ask spreads, returns, and compute trading costs
      df = submex_prefilt[ym == as.yearmon(names(x)[ind]),
                          .(cusip_id, ba, lba, ba_alt, lba_alt, retm)][df,,on = 'cusip_id']
      df[is.na(ba), ba := median(df[, ba], na.rm = T)]
      df[is.na(lba), lba := median(df[, lba], na.rm = T)]
      df[is.na(ba_alt), ba_alt := median(df[, ba_alt], na.rm = T)]
      df[is.na(lba_alt), lba_alt := median(df[, lba_alt], na.rm = T)]
      df[, .(weighted.mean(retm, w = w1, na.rm = T),
             sum(abs(w1-w0)),
             sum(ba/2*abs(w1-w0)),
             sum(ba_alt/2*abs(w1-w0)),
             nrow(df))]
      
    }
    
  })
  
  B = as.data.table(cbind(names(x)[seq], rbindlist(A)))
  names(B) = c('month', 'ret', 'trnvr', 'tcost', 'tcost_alt', 'nbonds')
  return(B)
  
}

# compute trading costs for all combinations of fyear and tc
tcost.full = lapply(comp, function(z) {tcost(x = z)})
names(tcost.full) = names(comp)

# compute net returns
res = lapply(tcost.full, function(x) {x[, c('netret', 'netret_alt') := .(ret - tcost, ret - tcost_alt)]})

# make xts's with results
REVxts = do.call(merge, lapply(res, function(x) {
  as.xts(x[, ret], order.by = as.yearmon(x[, month]))
  }))

REVxts_net = do.call(merge, lapply(res, function(x) {
  as.xts(x[, netret], order.by = as.yearmon(x[, month]))
  }))

REVxts_net_alt = do.call(merge, lapply(res, function(x) {
  as.xts(x[, netret_alt], order.by = as.yearmon(x[, month]))
}))

REVxtscum = do.call(merge, lapply(1+REVxts, cumprod))
REVxtscum_net = do.call(merge, lapply(1+REVxts_net, cumprod))
REVxtscum_net_alt = do.call(merge, lapply(1+REVxts_net_alt, cumprod))

REVxts_exmkt = do.call(merge, lapply(REVxts, function(x) {x - REVxts$mkt}))
REVxts_exmkt$mkt = NA

REVxts_net_exmkt = do.call(merge, lapply(REVxts_net, function(x) {x - REVxts_net$mkt}))
REVxts_net_exmkt$mkt = NA

REVxts_net_alt_exmkt = do.call(merge, lapply(REVxts_net_alt, function(x) {x - REVxts_net_alt$mkt}))
REVxts_net_alt_exmkt$mkt = NA

REVxts_exrf = do.call(merge, lapply(REVxts, function(x) {x - rf[index(x)]}))
REVxts_net_exrf = do.call(merge, lapply(REVxts_net, function(x) {x - rf[index(x)]}))
REVxts_net_alt_exrf = do.call(merge, lapply(REVxts_net_alt, function(x) {x - rf[index(x)]}))

# make a performance table

tab = cbind(12*colMeans(REVxts*100),
            sqrt(12)*sapply(REVxts*100, sd),
            sqrt(12)*colMeans(REVxts_exrf*100)/(sapply(REVxts_exrf*100, sd)),
            sqrt(12)*colMeans(REVxts_exmkt*100)/(sapply(REVxts_exmkt*100, sd)),
            12*colMeans(REVxts_net*100),
            sqrt(12)*sapply(REVxts_net*100, sd),
            sqrt(12)*colMeans(REVxts_net_exrf*100)/(sapply(REVxts_net_exrf*100, sd)),
            sqrt(12)*colMeans(REVxts_net_exmkt*100)/(sapply(REVxts_net_exmkt*100, sd)),
            12*colMeans(REVxts_net_alt*100),
            sqrt(12)*sapply(REVxts_net_alt*100, sd),
            sqrt(12)*colMeans(REVxts_net_alt_exrf*100)/(sapply(REVxts_net_alt_exrf*100, sd)),
            sqrt(12)*colMeans(REVxts_net_alt_exmkt*100)/(sapply(REVxts_net_alt_exmkt*100, sd)))
tab <- tab[-nrow(tab), ] # remove the market from the printout
colnames(tab) = rep(c('Mean', 'S.D.', 'SR', 'IR'), 3)
rownames(tab) = c(c('Baseline',  '~~Many funds', '~~Few funds',
                        '~~Big issuers', '~~Small issuers', '~~Inv. grade', '~~High yield'),
                  c('Baseline~',  '~~Many funds~', '~~Few funds~',
                    '~~Big issuers~', '~~Small issuers~', '~~Inv. grade~', '~~High yield~'))


## print summary pricing table for slides
addtorow = list()
addtorow$pos = list(-1, 0, 7)
addtorow$command = c("\\hline \\hline \\multicolumn{1}{l|}{} & \\multicolumn{4}{c|}{Before T-cost} &
                              \\multicolumn{4}{c|}{After T-cost (avg. bid-ask)} &
                              \\multicolumn{4}{c}{After T-cost (roundtrip)} \\\\\n",
                     "\\hline \\multicolumn{13}{c}{(A) Reversal: univariate sort on month $t-1$ return} \\\\\n",
                     "\\hline \\multicolumn{13}{c}{(B) Reversal: bi-variate sort on month $t-1$ return and rating} \\\\\n")

print(xtable(tab, auto = F,
             align = c('l|', rep('c',3), 'c|', rep('c',3), 'c|', rep('c', 4)),
             digits = rep(2,13)),
      file = '../Apps/ShareLaTeX/Bond reversals/performance_tab.tex', floating = F, include.colnames = T,
      include.rownames = T, size = 'footnotesize', sanitize.text.function=function(x){x},
      add.to.row = addtorow, hline.after = c(nrow(tab), nrow(tab)))


# plot stuff
#par(mfrow = c(1,2))
#plot.xts(REVxtscum, legend.loc = 'topleft')
#plot.xts(REVxtscum_net, legend.loc = 'topleft')

# plot stuff

#namef = '../Apps/ShareLaTeX/Bond reversals/invest_funds_gross.pdf'

#pdf(file = namef, family = 'CM Roman', width = 8, height = 8, pointsize = 24, onefile = F)

#plot(REVxtscum[, c('total', 'fewfunds', 'manyfunds', 'mkt')], major.ticks = 'years',
#     minor.ticks = NULL, main = 'Gross',
#     yaxis.right = FALSE, grid.ticks.on = 'years', grid.col = 'lightgray', grid.ticks.lty = 2,
#     lwd = 2, format.labels = '%b\n%Y', lty  = 1:4, col = 'black', on = NA)

#addLegend("topleft", legend.names = c('Long reversal (LR)', 'LR: few funds', 'LR: many funds',
#                                      'Market'), 
#          lty = 1:4, lwd= rep(3, 4), col = rep('black', 4), on = 0)

#dev.off()
#embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))


#namef = '../Apps/ShareLaTeX/Bond reversals/invest_funds_net.pdf'

#pdf(file = namef, family = 'CM Roman', width = 8, height = 8, pointsize = 24, onefile = F)

#plot(REVxtscum_net[, c('total', 'fewfunds', 'manyfunds', 'mkt')], major.ticks = 'years',
#     minor.ticks = NULL, main = 'Net of trading costs',
#     yaxis.right = FALSE, grid.ticks.on = 'years', grid.col = 'lightgray', grid.ticks.lty = 2,
#     lwd = 2, format.labels = '%b\n%Y', lty  = 1:4, col = 'black', on = NA)

#addLegend("bottomright", legend.names = c('Long reversal (LR)', 'LR: few funds', 'LR: many funds',
#                                          'Market'), 
#          lty = 1:4, lwd= rep(3, 4), col = rep('black', 4), on = 0)

#dev.off()
#embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))
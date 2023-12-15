# ==================================================================================================
# [] PREPARE DATA FOR REVERSAL ESTMATION
#  - RUN TIME: 15 mins
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer); library(lfe); library(plm); library(extrafont)
library(sandwich); library(lmtest); library(car); library(plotrix)

load_rmetrics_calendars(seq(2004,2020,1))

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')
source('./codes/ratfunc.R')
source('./codes/corstars.R')

# --------------------------------------------------------------------------------------------------
# (1) ADD VARIABLES
# --------------------------------------------------------------------------------------------------

# load data
load(file = './data_input/subtrace.RData')

# add coupon rate
c_rate = submergent[, .(cusip_id, COUPON)]
names(c_rate)[2] = 'c_rate'
subtrace = c_rate[subtrace,, on = 'cusip_id']
remove(c_rate)

# add call feature
call = submergent[, .(cusip_id, as.numeric(REDEEMABLE == 'Y'))]
names(call)[2] = 'call'
subtrace = call[subtrace,, on = 'cusip_id']
remove(call)

# make year-month and year-quarter variables 
subtrace[, ym := as.yearmon(trd_exctn_dt)]
subtrace[, yq := as.yearqtr(trd_exctn_dt)]
subtrace[, yr := year(trd_exctn_dt)]

# add dealers
load('./data_input/ndealers.RData')
names(ndealers)[3] = 'dealers'
ndealers[, ym := as.yearmon(month)]
ndealers[, month := NULL]
subtrace = ndealers[, .(cusip_id, ym, dealers)][subtrace,,on = c('cusip_id', 'ym')]
remove(ndealers)

# add number of funds
load('./data_input/nfunds.RData')
subtrace = nfunds[subtrace,, on = c('cusip_id', 'trd_exctn_dt')]
remove(nfunds)
gc()

## add equity
load(file = './data_input/subtrace_daily_traded.RData')
eqdf = subtrace_traded
remove(subtrace_traded)
eqdf[, c('issue_id', 'issuer') := NULL]
subtrace = eqdf[subtrace,, on = c('cusip_id', 'trd_exctn_dt')]
remove(eqdf)
gc()

# scale stuff, remove what is not needed
subtrace[, VE := VE/10^6]

# --------------------------------------------------------------------------------------------------
# (2) COMPUTE RETURNS, APPLY FILTERS
# --------------------------------------------------------------------------------------------------

setkey(subtrace, cusip_id, trd_exctn_dt)

# remove obs without volume
subtrace <- delete(subtrace, which(!is.na(subtrace$ivol)))
gc()

## compute returns
setkey(subtrace, cusip_id, trd_exctn_dt)
subtrace[, ret := (vwap+ai+coupon)/(shift(vwap) + shift(ai))-1, cusip_id]
subtrace[, plogret := log(vwap) - shift(log(vwap)), cusip_id]

submex_init = subtrace[, .(mean(size_hist, na.rm = T),
                           mean(rating_num, na.rm = T),
                           mean(nfunds, na.rm = T),
                           mean(cvol, na.rm = T),
                           mean(avg_ba, na.rm = T),
                           prod(1+ret)-1),
                       .(cusip_id, ym)]
names(submex_init) = c('cusip_id', 'ym', 'size',
                       'rat', 'nfunds', 'cvol', 'avg_ba', 'retm')

## apply filters
subtrace = subtrace[rating_num <= 21 &
                      age > 0     &
                      spread > 0 &
                      size_hist > 0    &
                      mat > 1     &
                      trd_exctn_dt >= as.IDate('2005-01-01') &
                      trd_exctn_dt <= as.IDate('2018-12-31') &
                      !is.na(dealers)]

# --------------------------------------------------------------------------------------------------
# (3) WINSORIZE
# --------------------------------------------------------------------------------------------------

# winsorize volume
volvars = c('vol', 'volb', 'vols', 'vold', 'volexd', 'cvol',
            'cumvol', 'cumvolb', 'cumvols', 'cumvold')
subtrace[, c('vol', 'volb', 'vols', 'vold', 'volexd', 'cvol',
             'cumvol', 'cumvolb', 'cumvols', 'cumvold') := lapply(.SD, function(x) {
               winsorize(x, 0, 0.99)
             }),,
         .SDcols = volvars]
subtrace[, ivol := winsorize(ivol, 0.01, 0.99)]

# winsorize spread, avg ba, and total volume
subtrace[, spread := winsorize(spread, 0, 0.999)]
subtrace[, avg_ba := winsorize(avg_ba, 0, 0.999)]

## winsorize returns
subtrace[, ret := winsorize(ret, 0.001, 0.999)]
subtrace[, plogret := winsorize(plogret, 0.001, 0.999)]

## compute days since last trade
subtrace[, pdate := shift(trd_exctn_dt), cusip_id]
subtrace[!is.na(pdate), dsince := bizdays(pdate, trd_exctn_dt, cal = 'Rmetrics/NYSE'), cusip_id]
subtrace[, pdate := NULL]

## remove observations with missing dates since last trade (first obs per bond)
subtrace = subtrace[!is.na(dsince)]

## remove bonds with a single observation
subtrace = subtrace[!cusip_id %in% subtrace[, .(.N), cusip_id][N==1, cusip_id]]

# --------------------------------------------------------------------------------------------------
# (4) COMPUTE ADDITIONAL VOLUME METRICS
# --------------------------------------------------------------------------------------------------

# compute volume measures
wdth = 60

subtrace[, cvolmean := (if (.N < wdth) {
  rep(-9999, .N)
} else {
  rollapply(cvol, mean, width = wdth, fill = NA, align = 'right', na.rm = T)
}), cusip_id]

subtrace[, absivolmean := (if (.N < wdth) {
  rep(-9999, .N)
} else {
  rollapply(abs(ivol), mean, width = wdth, fill = NA, align = 'right', na.rm = T)
}), cusip_id]

subtrace[cvolmean == -9999, cvolmean := NA]
subtrace[absivolmean == -9999, absivolmean := NA]

subtrace[, ctrnvr := log(cvol + 0.00000255) - log(cvolmean + 0.00000255)]
subtrace[, ctrnvr_alt := log((cvol+0.00000255)^2)]

subtrace[, itrnvr := log(abs(ivol) + 0.00000255) - log(absivolmean + 0.00000255)]
subtrace[, itrnvr_alt := log((abs(ivol)+0.00000255)^2)]

# --------------------------------------------------------------------------------------------------
# (5) SELECT SMOOTH PERIODS
# --------------------------------------------------------------------------------------------------

## index periods
indperiods = function(x, th) {
  
  ind = vector('numeric', length(x))
  ind[1] = 1
  
  if (length(x) > 1) {
    
    for (i in 2:length(x)) {
      if ( ((x[i] <= th) & (x[i-1] <= th)) | ((x[i] <= th) & is.na(x[i-1])) ) {
        ind[i] = ind[i-1] 
      } else {
        ind[i] = ind[i-1] + 1
      }
    }
    
  }
  
  return(ind)
  
}

## apply indices daily, for each cusip and year-quarter
subtrace[, ind := indperiods(x = dsince, th = 3), .(cusip_id, yq)]

## compute number of observations per period per bond per year-quarter
subtrace[, Nind := .N, .(cusip_id, yq, ind)]

## leave only chunks with more than Nlim observations
subtrace_prefilt = subtrace
Nlim = 40
subtrace = subtrace_prefilt[Nind >= Nlim]

## remove crossing of IG-HY boundary
subtrace[, cross := (sum(rating_num <= 10)>0 & sum(rating_num > 10)>0),  .(cusip_id, yq, ind)]
subtrace = subtrace[cross == F]

## compute additional volume indicators
subtrace[, ctrnvr_alt_2 :=  cvol - mean(cvol), .(cusip_id, yq, ind)]
subtrace[, itrnvr_alt_2 :=  ivol/sd(ivol), .(cusip_id, yq, ind)]
subtrace[, ctrnvr_alt_3 :=  (cvol - mean(cvol))/sd(cvol), .(cusip_id, yq, ind)]
subtrace[, itrnvr_alt_3 :=  (abs(ivol) - mean(abs(ivol)))/sd(abs(ivol)), .(cusip_id, yq, ind)]
subtrace[, itrnvr_alt_4 :=  (ivol - mean(ivol))/sd(ivol), .(cusip_id, yq, ind)]

## compute signed volumes and price impacts of the order flow
subtrace[, slogivol := sign(ivol)*log(abs(ivol*size_hist*10^6)+1)]

## winsorize additional volume indicators
subtrace[, ctrnvr_alt_3 := winsorize(ctrnvr_alt_3, 0, 0.99)]
subtrace[, itrnvr_alt_3 := winsorize(itrnvr_alt_3, 0, 0.99)]
subtrace[, itrnvr_alt_4 := winsorize(itrnvr_alt_4, 0, 0.99)]

## compute return volatility
subtrace[, sigma := sd(ret, na.rm = T), .(cusip_id, yq, ind)]

# median length of the smooth period
unique(subtrace[, .(cusip_id, yq, ind, Nind)])[, median(Nind)]

# remove observations without turnover
subtrace = subtrace[!is.na(ctrnvr_alt_3)]

# --------------------------------------------------------------------------------------------------
# (6) ADD EQUITY DATA PER SMOOTH PERIOD AND SAVE EVERYTHING
# --------------------------------------------------------------------------------------------------

setkey(subtrace, cusip_id, trd_exctn_dt)

# calculate equity volatility within each active period
subtrace[, SE := sd(RET, na.rm = T), .(cusip_id, yq, ind)]

save(subtrace, subtrace_prefilt, submex_init, submergent, Nlim,
     file = './data_input/subtrace_with_smoothper_peryq.RData')

# --------------------------------------------------------------------------------------------------
# (7) PRINT SUMMARY STATS
# --------------------------------------------------------------------------------------------------

# run some summary stats
sumtabfunc = function(x) {
  
  tab = rbind(descfull(x[, size_hist], dig = 0),
              descfull(x[, mat]),
              descfull(x[, c_rate]),
              descfull(x[, rating_num]),
              descfull(x[, age]),
              descfull(x[, ret*100]),
              descfull(x[, cvol]),
              descfull(x[, ivol]),
              descfull(x[, abs(ivol)]),
              descfull(x[, avg_ba]),
              descfull(x[, nfunds], dig = 1),
              descfull(x[, dealers], dig = 1),
              descfull(x[, VE], dig = 1),
              descfull(x[, stock_ba]),
              descfull(x[, distea], dig = 1))
  
  rownames(tab) = c('Issue size, mln \\$', 'Maturity, years', 'Coupon rate, \\%',
                    'Rating', 'Age, years', 'Total bond return, \\%', 
                    'C-to-C volume, \\% of size', 'C-to-D volume, \\% of size',
                    '$|\\text{C-to-D volume}|$, \\% of size', 'Realized bond bid-ask, \\%', 
                    'No. mutual fund owners', 'No. dealers', 'Issuer equity value, bln \\$',
                    'Stock bid-ask, \\%', 'Days to earnings announcement'
  )
  tab
  
}

tab = sumtabfunc(subtrace)

print(xtable(tab, align = c('l|', rep('c', 2), 'c|', rep('c', 6), '|c'),
             auto = F), file = '../Apps/ShareLaTeX/Bond reversals/sumstat_peryq.tex', floating = F,
      include.colnames = T, include.rownames = T, size = 'footnotesize',
      sanitize.text.function=function(x){x}, hline.after = c(-1, 0, nrow(tab)))

remove(tab)
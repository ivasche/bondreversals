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

# make year-month variable 
subtrace[, ym := as.yearmon(trd_exctn_dt)]

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

# add alt prices
load(file = './data_input/subprices.RData')
subprices = unique(subprices[, .(cusip_id, trd_exctn_dt, vwmp, vwapb, vwapb_alt)])
gc()

# clean
subprices[, Nobs := .N, .(cusip_id, trd_exctn_dt)]
subprices[, tokeep := F]
subprices[Nobs == 1, tokeep := T]
subprices[Nobs > 1 & is.na(vwapb) == F & is.na(vwapb_alt) == F, tokeep := T]
subprices <- delete(subprices, subprices$tokeep)

# merge
subtrace <- subprices[subtrace,, on = c('cusip_id', 'trd_exctn_dt')]
remove(subprices)
gc()

## compute returns
setkey(subtrace, cusip_id, trd_exctn_dt)
subtrace[, ret := (vwap+ai+coupon)/(shift(vwap) + shift(ai))-1, cusip_id]
subtrace[, plogret := log(vwap) - shift(log(vwap)), cusip_id]

subtrace[, ret_mid := (vwmp+ai+coupon)/(shift(vwmp) + shift(ai))-1, cusip_id]
subtrace[, plogret_mid := log(vwmp) - shift(log(vwmp)), cusip_id]

subtrace[, ret_big := (vwapb+ai+coupon)/(shift(vwapb) + shift(ai))-1, cusip_id]
subtrace[, plogret_big := log(vwapb) - shift(log(vwapb)), cusip_id]

subtrace[, ret_big_alt := (vwapb_alt+ai+coupon)/(shift(vwapb_alt) + shift(ai))-1, cusip_id]
subtrace[, plogret_big_alt := log(vwapb_alt) - shift(log(vwapb_alt)), cusip_id]

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
gc()

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

subtrace[, ret_mid := winsorize(ret_mid, 0.001, 0.999)]
subtrace[, plogret_mid := winsorize(plogret_mid, 0.001, 0.999)]

subtrace[, ret_big := winsorize(ret_big, 0.001, 0.999)]
subtrace[, plogret_big := winsorize(plogret_big, 0.001, 0.999)]

subtrace[, ret_big_alt := winsorize(ret_big_alt, 0.001, 0.999)]
subtrace[, plogret_big_alt := winsorize(plogret_big_alt, 0.001, 0.999)]

## compute days since last trade
subtrace[, pdateb := shift(trd_exctn_dt), cusip_id]
subtrace[, flagb := shift(vwapb), cusip_id]
subtrace[is.na(flagb), pdateb := NA]
subtrace[!is.na(pdateb) & !is.na(vwapb),
         dsince := bizdays(pdateb, trd_exctn_dt, cal = 'Rmetrics/NYSE'), cusip_id]
subtrace[, c('pdateb', 'flagb') := NULL]

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
  
  for (i in 2:length(x)) {
    if ( ((x[i] <= th) & (x[i-1] <= th)) | ((x[i] <= th) & is.na(x[i-1])) ) {
      ind[i] = ind[i-1] 
    } else {
      ind[i] = ind[i-1] + 1
    }
  }
  ind
  
}

## apply indices daily
setkey(subtrace, cusip_id, trd_exctn_dt)
subtrace[, ind := indperiods(x = dsince, th = 3), cusip_id]
#subtrace[, ind_alt := indperiods(x = dsince_alt, th = 3), cusip_id]

## compute number of observations per period per bond
subtrace[, Nind := .N, .(cusip_id, ind)]

## leave only chunks with more than Nlim observations
subtrace_prefilt = subtrace
Nlim = 60
subtrace = subtrace_prefilt[Nind >= Nlim]
gc()

## remove crossing of IG-HY boundary
subtrace[, cross := (sum(rating_num <= 10)>0 & sum(rating_num > 10)>0),  .(cusip_id, ind)]
subtrace <- delete(subtrace, keep.idxs = (subtrace$cross == F))
gc()

## compute additional volume indicators
subtrace[, ctrnvr_alt_2 :=  cvol - mean(cvol), .(cusip_id, ind)]
subtrace[, itrnvr_alt_2 :=  ivol/sd(ivol), .(cusip_id, ind)]
subtrace[, ctrnvr_alt_3 :=  (cvol - mean(cvol))/sd(cvol), .(cusip_id, ind)]
subtrace[, itrnvr_alt_3 :=  (abs(ivol) - mean(abs(ivol)))/sd(abs(ivol)), .(cusip_id, ind)]
subtrace[, itrnvr_alt_4 :=  (ivol - mean(ivol))/sd(ivol), .(cusip_id, ind)]

## compute signed volumes and price impacts of the order flow
subtrace[, slogivol := sign(ivol)*log(abs(ivol*size_hist*10^6)+1)]

## winsorize additional volume indicators
subtrace[, ctrnvr_alt_3 := winsorize(ctrnvr_alt_3, 0, 0.99)]
subtrace[, itrnvr_alt_3 := winsorize(itrnvr_alt_3, 0, 0.99)]
subtrace[, itrnvr_alt_4 := winsorize(itrnvr_alt_4, 0, 0.99)]

## yet another volume test
subtrace[, ctrnvr_alt_5 :=  (log(cvol + 1)-mean(log(cvol + 1)))/sd(log(cvol + 1)), .(cusip_id, ind)]
subtrace[, itrnvr_alt_5 :=  (log(abs(ivol) + 1)-mean(log(abs(ivol) + 1)))/sd(log(abs(ivol) + 1)), .(cusip_id, ind)]

## winsorize additional volume indicators
subtrace[, ctrnvr_alt_5 := winsorize(ctrnvr_alt_5, 0, 0.99)]
subtrace[, itrnvr_alt_5 := winsorize(itrnvr_alt_5, 0, 0.99)]

## compute bond return volatility
subtrace[, sigma := sd(ret, na.rm = T), .(cusip_id, ind)]

# median length of the smooth period
unique(subtrace[, .(cusip_id, ind, Nind)])[, median(Nind)]
gc()

# remove observations without turnover
subtrace <- delete(subtrace, keep.idxs = !is.na(subtrace$ctrnvr_alt_3))
gc()

# --------------------------------------------------------------------------------------------------
# (6) ADD EQUITY DATA PER SMOOTH PERIOD AND SAVE EVERYTHING
# --------------------------------------------------------------------------------------------------

setkey(subtrace, cusip_id, trd_exctn_dt)

# calculate equity volatility within each active period
subtrace[, SE := sd(RET, na.rm = T), .(cusip_id, ind)]

subtrace_large <- subtrace
remove(subtrace)
gc()

save(subtrace_large, 
     file = './data_input/subtrace_with_smoothper_large.RData')

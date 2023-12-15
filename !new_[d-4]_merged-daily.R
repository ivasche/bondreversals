# ==================================================================================================
# [5] PREpARE DATASET FOR REVERSAL ESTMATION
#  - RUN TIME: 10 min
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer)

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')
source('./codes/ratfunc.R')

# --------------------------------------------------------------------------------------------------
# (1) LOAD DATA
# --------------------------------------------------------------------------------------------------

# load data
load('./data_input/subtrace_daily_outstanding_rating.RData')
setkey(subtrace, cusip_id, trd_exctn_dt)

# remove non-business days
load_rmetrics_calendars(seq(2004, 2020, 1))
bdays = bizseq(from = min(subtrace$trd_exctn_dt),
               to = max(subtrace$trd_exctn_dt), cal = 'Rmetrics/NYSE')
bdays = as.IDate(bdays)
subtrace <- delete(subtrace, keep.idxs = which(subtrace$trd_exctn_dt %fin% bdays))
gc()

# clean spreads
subtrace[ytm == -999900, ytm := NA]
subtrace[ytm_rf == -999900, ytm_rf := NA]
#subtrace[dur == -9999, dur := NA]

# compute age
subtrace[, age := as.yearmon(trd_exctn_dt) - as.yearmon(issue_date)]

# clean
subtrace[, c('cusip6', 'ISSUE_ID', 'issue_date', 'dated_date', 'first_int_date',
             'last_int_date', 'int_freq') := NULL]
gc()

# --------------------------------------------------------------------------------------------------
# (2) ADD LIQUIDITY
# --------------------------------------------------------------------------------------------------

# read liquidity data
load('./data_input/liquidity_daily.RData')
gc()

# add
subtrace = dliq[subtrace[cusip_id %in% unique(dliq[, cusip_id])],, on = c('cusip_id', 'trd_exctn_dt')]
remove(dliq)
gc()

# read volume data
load('./data_input/volume_daily.RData')
voldf[, yq := as.yearqtr(trd_exctn_dt)]
voldf[, td := (.N - sum(vol==0))/.N, .(cusip_id, yq)] # % of business days with trades, per quarter
voldf[, yq := NULL]
subtrace = voldf[subtrace,, on = c('cusip_id', 'trd_exctn_dt')]
remove(voldf)
gc()

## add roundtrip volume
load(file = './data_input/vlm_irtc.RData')
names(rtvlm)[3] <- 'rtv'
rtvlm[, trd_exctn_dt := as.IDate(as.Date(as.character(trd_exctn_dt), format = '%Y%m%d'))]
subtrace <- rtvlm[subtrace,, on = c('cusip_id', 'trd_exctn_dt')]
subtrace[is.na(rtv), rtv := 0]
remove(rtvlm)
gc()

setkey(subtrace, cusip_id, trd_exctn_dt)

# make sure that BUYs are client buys (= diller SELLs)
subtrace[, irow := .I]
subtrace[, `:=` ('vol_s' = vols, 'vol_b' = volb, 'vol_d'  = vold)]
subtrace[, volb   := vol_s]
subtrace[, vols   := vol_b]
subtrace[, vold   := vol_d]
subtrace[, volexd := volb + vols]
subtrace[, cvol   := min(volb, vols), irow]
subtrace[, irow := NULL]
subtrace[, ivol   := volb - vols]
subtrace[, c('vol_b', 'vol_d', 'vol_s') := NULL]

subtrace[, cumvols := b_cumvol]
subtrace[, cumvolb := s_cumvol]
subtrace[, cumvold := d_cumvol]
subtrace[, c('b_cumvol', 's_cumvol', 'd_cumvol') := NULL]
gc()

# rescale volumes
volvars = c('vol', 'volb', 'vols', 'vold', 'volexd', 'cvol', 'ivol', 'rtv',
            'cumvol', 'cumvolb', 'cumvols', 'cumvold')

subtrace[, c('vol', 'volb', 'vols', 'vold', 'volexd', 'cvol', 'ivol', 'rtv',
             'cumvol', 'cumvolb', 'cumvols', 'cumvold') := lapply(.SD, function(x) {x/10^3/size_hist}),,
         .SDcols = volvars]
gc()

# compute spreads and the share of dealer volume
subtrace[, spread := ytm - ytm_rf]
subtrace[, dshare := vold/vol*100]
subtrace[, dcumshare := cumvold/cumvol*100]

# cut submergent
submergent = submergent[cusip_id %fin% unique(subtrace$cusip_id)]

# save dataset
save(subtrace, submergent, file = './data_input/subtrace.RData')

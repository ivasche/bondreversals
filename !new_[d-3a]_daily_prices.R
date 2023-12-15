# ==================================================================================================
# [3a] COMPUTE ALL SORTS OF DAILY PRICES
#  - RUN TIME: 15 MINS
# ==================================================================================================

tick = Sys.time()

# load packages
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime)
library(xtable); library(stargazer)
library(haven)

# set locale
Sys.setlocale(locale = 'US')

# load user-defined function
source('./codes/delete.R')
source('./codes/ratfunc.R')
source('./codes/desc.R')

# --------------------------------------------------------------------------------------------------
# read trace records, filter execution time
# --------------------------------------------------------------------------------------------------

load('./data_input/subtrace_tick.RData')
trace = subtrace
remove(subtrace)
gc()

## mark small trades (first definition)
trace[, small100 := F]
trace[entrd_vol_qt <= 10000, small100 := T]

## mark small trades (second definition)
trace[, small100alt := F]
trace[entrd_vol_qt <= 100000, small100alt := T]

## compute vwap, vwmp (volume-wighted mid price) and vwap without small trades
trace[, vwap := weighted.mean(x = rptd_pr, w = entrd_vol_qt), .(cusip_id, trd_exctn_dt)]
gc()

trace[, vwp := weighted.mean(x = rptd_pr, w = entrd_vol_qt), .(cusip_id, trd_exctn_dt, rpt_side_cd)]
trace[, vwmp := mean(unique(vwp)), .(cusip_id, trd_exctn_dt)]
trace[, vwp := NULL]
gc()

trace[small100 == F, vwapb := weighted.mean(x = rptd_pr, w = entrd_vol_qt), .(cusip_id, trd_exctn_dt)]
trace[small100alt == F, vwapb_alt := weighted.mean(x = rptd_pr, w = entrd_vol_qt), .(cusip_id, trd_exctn_dt)]
gc()

trace[, c('company_symbol', 'entrd_vol_qt', 'rptd_pr', 'rpt_side_cd', 'small100', 'small100alt') := NULL]
gc()

subprices = unique(trace[, .(cusip_id, trd_exctn_dt, vwap, vwmp, vwapb, vwapb_alt)])
remove(trace)
subprices[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]
setkey(subprices, cusip_id, trd_exctn_dt)

save(subprices, file = './data_input/subprices.RData')
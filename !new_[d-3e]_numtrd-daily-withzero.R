# ==================================================================================================
# [3e] COMPUTE NUMBER OF NO-TRADING DAYS FOR DAILY DATA INCL. NO TRADING DAYS!
#  - RUN TIME: 
# ==================================================================================================

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
# CALCULATE NUMBER OF TRADES
# --------------------------------------------------------------------------------------------------

load('./data_input/subtrace_tick.RData')
trace = subtrace
remove(subtrace, submergent)
trace[, company_symbol := NULL]
gc()

## number of trades
numdf = trace[, .(.N), .(cusip_id, trd_exctn_dt, rpt_side_cd)]
numdf = dcast(numdf, cusip_id +  trd_exctn_dt ~ rpt_side_cd, value.var = 'N', fill = 0)
names(numdf)[3:5] = paste0('num', tolower(names(numdf)[3:5]))
gc()

# clean memory
remove(trace)
gc()

#convert date field
numdf[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]
gc()

# --------------------------------------------------------------------------------------------------
# MERGE
# --------------------------------------------------------------------------------------------------

# load back number of trades (incl missing)
load(file = './data_input/numtrd_interim.RData')
gc()

# merge 
numtrd = rbindlist(numtrd)
gc()

numtrd = numdf[numtrd,,on = c('cusip_id', 'trd_exctn_dt')]
remove(numdf)

# replace NAs (no trading dates) with zero volumes
numtrd[is.na(numb), numb := 0]
numtrd[is.na(numd), numd := 0]
numtrd[is.na(nums), nums := 0]

# calculate cumsum
numtrd[, c('b_cumzero', 's_cumzero', 'd_cumzero', 'cumzero') := lapply(.SD, function(x) {
  
  Reduce(`+`, shift(x==0, 0:30))
  
}), cusip_id, .SDcols = c('numb', 'nums', 'numd', 'trds')]

numdf <- numtrd
remove(numtrd)
gc()

save(numdf, file = paste0('./data_input/numtrd_daily.RData'))
file.remove('./data_input/numtrd_interim.RData')
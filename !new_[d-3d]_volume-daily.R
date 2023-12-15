# ==================================================================================================
# [3d] COMPUTE VOLUME INDICATORS FOR DAILY DATA INCL. NO TRADING DAYS!
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
# CALCULATE VOLUME
# --------------------------------------------------------------------------------------------------

load('./data_input/subtrace_tick.RData')
trace = subtrace
remove(subtrace, submergent)
trace[, company_symbol := NULL]
gc()

## VOLUME
voldf = trace[, sum(entrd_vol_qt), .(cusip_id, trd_exctn_dt, rpt_side_cd)]
voldf = dcast(voldf, cusip_id + trd_exctn_dt ~ rpt_side_cd, value.var = 'V1', fill = 0)
names(voldf)[3:5] = paste0('vol', tolower(names(voldf)[3:5]))
gc()

# clean memory
remove(trace)
gc()

#convert date field
voldf[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]
gc()

# winsorize

# --------------------------------------------------------------------------------------------------
# MERGE
# --------------------------------------------------------------------------------------------------

# load back number of trades (incl missing)
load(file = './data_input/numtrd_interim.RData')
gc()

# merge 
numtrd = rbindlist(numtrd)
gc()

numtrd = voldf[numtrd,,on = c('cusip_id', 'trd_exctn_dt')]
remove(voldf)

# replace NAs (no trading dates) with zero volumes
numtrd[is.na(volb), volb := 0]
numtrd[is.na(vold), vold := 0]
numtrd[is.na(vols), vols := 0]
numtrd[, vol := volb+vold+vols]

# calculate cumsum
numtrd[, c('b_cumvol', 's_cumvol', 'd_cumvol', 'cumvol') := lapply(.SD, function(x) {
  
  Reduce(`+`, shift(x, 0:30))
  
  }), cusip_id, .SDcols = c('volb', 'vols', 'vold', 'vol')]

numtrd[, trds := NULL]
voldf <- numtrd
remove(numtrd)
gc()

save(voldf, file = paste0('./data_input/volume_daily.RData'))
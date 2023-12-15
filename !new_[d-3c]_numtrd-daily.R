# ==================================================================================================
# [3c] COMPUTE NO TRADING DAYS
#  - RUN TIME: 15 mins
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
# COMPUTE NUMBER OF TRADES WITH MISSING DATES
# --------------------------------------------------------------------------------------------------

load('./data_input/subtrace_tick.RData')
trace = subtrace
remove(subtrace)
trace[, company_symbol := NULL]
gc()

# save date boundaries
firstdate = as.IDate(as.character(min(trace[, trd_exctn_dt])), format = '%Y%m%d')
lastdate = as.IDate(as.character(max(trace[, trd_exctn_dt])), format = '%Y%m%d')

## record the number of trades per actual trading day
toadd = trace[, length(rptd_pr), .(cusip_id, trd_exctn_dt)]
names(toadd)[3] = 'trds'
toadd[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]
remove(trace)
gc()

# leave in submergent only the necessary stuff
datecols = c("MATURITY", "OFFERING_DATE")
submergent[, colnames(submergent)[!colnames(submergent) %fin% c('cusip_id', datecols)] := NULL]
submergent[, OFFERING_DATE := as.IDate(as.character(OFFERING_DATE), format = '%Y%m%d')]
submergent[, MATURITY := as.IDate(as.character(MATURITY), format = '%Y%m%d')]

## add info on issuance date and maturity
numtrd  = submergent
remove(submergent)
toremove = c(numtrd[MATURITY < firstdate][, cusip_id], numtrd[OFFERING_DATE > lastdate, cusip_id], 
             numtrd[is.na(OFFERING_DATE) | is.na(MATURITY), cusip_id],
             numtrd[OFFERING_DATE >= MATURITY, cusip_id])
numtrd = numtrd[!cusip_id %fin% toremove]
gc()

## construct a sequence of dates for each bond between [max(issuance, sample start), min(maturity, sample end)]
numtrd  = split(numtrd , numtrd$cusip_id)

## create a sequence of business days
load_rmetrics_calendars(seq(2004, 2020, 1))
bdays = bizseq(from = firstdate, to = lastdate, cal = 'Rmetrics/NYSE')
bdays = as.IDate(bdays)

numtrd  = pblapply(numtrd , function(x) {
  
  A = data.table(trd_exctn_dt = seq(max(x$OFFERING_DATE, firstdate), min(x$MATURITY, lastdate), by = "1 day"))
  return(A[trd_exctn_dt %fin% bdays])
  
  })
numtrd  = rbindlist(numtrd , idcol =  'cusip_id')
gc()

# add info on which dates there're trades
numtrd  = toadd[numtrd,, on = c('cusip_id', 'trd_exctn_dt')]
remove(toadd)
gc()

# clean
numtrd[is.na(trds), trds := 0]
numtrd[cusip_id %fin% toremove, trds := NA]
gc()

numtrd = split(numtrd, numtrd$cusip_id)
save(numtrd, file = paste0('./data_input/numtrd_interim.RData'))
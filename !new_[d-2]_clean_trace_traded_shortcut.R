# ==================================================================================================
# [2] MERGE TRACE WITH COMPUSTAT AND CRSP (ONLY TRADED COMPANIES' DEBT)
#  - RUN TIME: ABOUT 40 MIN
# ==================================================================================================

# load packages
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime)
library(xtable); library(stargazer)
library(haven); library(nleqslv)

# set locale
Sys.setlocale(locale = 'US')

# load user-defined function
source('./codes/delete.R')
source('./codes/ratfunc.R')
source('./codes/desc.R')

# load a trace-crsp link
trace.crsp = fread('./data_input/trace_crsp_link.csv')
trace.crsp[, trace_startdt := as.IDate(as.character(trace_startdt), format = '%Y%m%d')]
trace.crsp[, trace_enddt := as.IDate(as.character(trace_enddt), format = '%Y%m%d')]
trace.crsp[, crsp_startdt := as.IDate(as.character(crsp_startdt), format = '%Y%m%d')]
trace.crsp[, crsp_enddt := as.IDate(as.character(crsp_enddt), format = '%Y%m%d')]
trace.crsp[, link_startdt := as.IDate(as.character(link_startdt), format = '%Y%m%d')]
trace.crsp[, link_enddt := as.IDate(as.character(link_enddt), format = '%Y%m%d')]

names(trace.crsp)[1] = 'cusip_id'

# load daily trace data
load('./data_input/subtrace_daily_outstanding_rating.RData')

# merge (will yield several permnos for some bonds on some days)
subtrace = trace.crsp[subtrace,, on = c('cusip_id'), allow.cartesian = T]
subtrace[trd_exctn_dt >= link_startdt & trd_exctn_dt <= link_enddt, ind := 1]
subtrace = delete(subtrace, keep.idxs = which(subtrace$ind == 1))

# break ties (the latest trace entry is a priority, else the longest trace sample)
breakties = function(traceend, tracestart) {
  
  ind <- (traceend == max(traceend, na.rm = T))
  if (sum(ind) > 1) {
        ind <- (tracestart == min(tracestart, na.rm = T))
  }
  return(ind)
  
}

subtrace[, ind := NULL]
subtrace[, ind := breakties(trace_enddt, trace_startdt), .(cusip_id, trd_exctn_dt)]
subtrace = delete(subtrace, which(subtrace$ind == 1))
gc()

subtrace[, c('trace_startdt', 'trace_enddt', 'crsp_startdt', 'crsp_enddt', 
             'link_startdt', 'link_enddt', 'ind') := NULL]
remove(trace.crsp)
gc()

# --------------------------------------------------------------------------------------------------
# MERGE TRACE WITH CRSP-based DATA AND COMPUSTAT FOR FURTHER DD CALCULATONS
# --------------------------------------------------------------------------------------------------

# read crsp daily data
crsp = fread('./data_input/crsp_only_new.csv')

# remove weird pricing, convert dates into date format and returns into numeric format
crsp = crsp[PRC > 0]
crsp[, date := as.IDate(as.character(date), format = '%Y%m%d')]
crsp[, RET := as.numeric(RET)]
crsp = crsp[!is.na(RET)]

# leave only permn(c)os from trace
crsp = crsp[PERMNO %fin% unique(subtrace[, PERMNO]) | PERMCO %fin% unique(subtrace[, PERMCO])]
names(crsp)[2] = 'trd_exctn_dt'

# compute VE = firm equity value
crsp[, VE := PRC*SHROUT]

# compute stock bid-ask
crsp[, stock_ba := (ASK-BID)/PRC*100]
crsp[, stock_ba := winsorize(stock_ba, 0.01, 0.99)]

# compute SE = returns volatility (equity volatility)
setkey(crsp, PERMNO, trd_exctn_dt)
crsp[, yq := as.yearqtr(trd_exctn_dt)]
crsp[, yr := year(trd_exctn_dt)]
#SE = function(x) {rollapply(as.xts(x$RET, order.by = as.Date(x$trd_exctn_dt)), FUN = sd,
#                            width = min(30, nrow(x)), na.rm = T)}
#crsp[, SE := SE(.SD), by = PERMNO]
crsp[, SEpq := sd(RET, na.rm = T), .(PERMNO, yq)]
crsp[, SEpy := sd(RET, na.rm = T), .(PERMNO, yr)]

# winsorize
crsp[, SEpq := winsorize(SEpq, lower = 0.01, upper = 0.99)]
crsp[, SEpy := winsorize(SEpy, lower = 0.01, upper = 0.99)]

# --------------------------------------------------------------------------------------------------
# MERGE WITH IBES
# --------------------------------------------------------------------------------------------------

# read ibes data and add cusips to surprise history
ibes_sue = fread('./data_input/ibes_sue.csv')
ibes_actuals = fread('./data_input/ibes_actuals.csv')

cus <- unique(ibes_actuals[, .(CUSIP, OFTIC, ANNDATS)])
names(cus)[3] = 'anndats'
ibes_sue <- cus[ibes_sue,,on = c('OFTIC', 'anndats')]
ibes_sue <- ibes_sue[!is.na(CUSIP)]
remove(cus, ibes_actuals)
gc()

# leave only securities from TRACE
ibes_sue <- ibes_sue[CUSIP %fin% c(unique(crsp[, CUSIP]), unique(crsp[, NCUSIP]))]

ibes_sue[, anndats := as.IDate(as.character(anndats), format = '%Y%m%d')]
ibes_sue[, trd_exctn_dt := anndats]

# create year-month var for the last month of a reporting quarter
# and for multiple reporting periods for an announcement date pick the most recent one
ibes_sue[, repmon := as.yearmon(paste0(PYEAR, '-', PMON))]
ibes_sue[, annmon := as.yearmon(anndats)]
ibes_sue[, diffrep := 12*(annmon - repmon)]
ibes_sue <- ibes_sue[diffrep > 0]
ibes_sue[, ind := (diffrep == diffrep[which.min(diffrep)]), .(CUSIP, anndats)]
ibes_sue <- delete(ibes_sue, ibes_sue$ind)
ibes_sue[, ind := NULL]

# add announcement days to crsp
crsp <- ibes_sue[, .(CUSIP, anndats, trd_exctn_dt)][crsp,, on = c('CUSIP', 'trd_exctn_dt'), roll = -Inf]

alt <- ibes_sue[, .(CUSIP, anndats, trd_exctn_dt)]
names(alt)[1:2] = c('NCUSIP', 'anndats_alt')

crsp <- alt[crsp,, on = c('NCUSIP', 'trd_exctn_dt'), roll = -Inf]
crsp[is.na(anndats) & !is.na(anndats_alt), anndats := anndats_alt]
crsp[, anndats_alt := NULL]
remove(alt)
gc()

# mark earnings announcement days
crsp[!is.na(anndats), eaday := 1]
crsp[is.na(eaday), eaday := 0]

# distance to earnings announcment
crsp[, distea := anndats - trd_exctn_dt]
crsp[distea > 92, c('distea', 'anndats', 'eaday') := NA]

# add analyst dispersion
crsp <- ibes_sue[, .(CUSIP, trd_exctn_dt, suescore, surpstdev)][crsp, on = c('CUSIP', 'trd_exctn_dt'), roll = +Inf]
alt <- ibes_sue[, .(CUSIP, trd_exctn_dt, suescore, surpstdev)]
names(alt) = c('NCUSIP', 'trd_exctn_dt', 'suescore_alt', 'surpstdev_alt')
crsp <- alt[crsp,, on = c('NCUSIP', 'trd_exctn_dt'), roll = +Inf]
crsp[is.na(suescore) & !is.na(suescore_alt), suescore := suescore_alt]
crsp[is.na(surpstdev) & !is.na(surpstdev_alt), surpstdev := surpstdev_alt]
crsp[, c('suescore_alt', 'surpstdev_alt') := NULL]
remove(alt)
gc()

# scale surpstdev
crsp[, surpstdev := surpstdev/PRC*100]

# winsorize
crsp[, suescore := winsorize(suescore, 0.01, 0.99)]
crsp[, surpstdev := winsorize(surpstdev, 0, 0.99)]

# merge
subtrace = crsp[, .(PERMNO, trd_exctn_dt, PRC, RET, VOL, stock_ba, eaday, distea, suescore, surpstdev,
                    SHROUT, VE, SEpq, SEpy)][subtrace,, on = c('PERMNO', 'trd_exctn_dt')]

subtrace_traded = subtrace[, .(cusip_id, trd_exctn_dt, issue_id, issuer, PERMCO, PERMNO, 
                               PRC, RET, VOL, stock_ba, eaday, distea, suescore, surpstdev, SHROUT, VE, SEpq, SEpy)]
subtrace_traded <- delete(subtrace_traded, !is.na(subtrace_traded$RET))
remove(subtrace)
gc()


save(subtrace_traded, file = './data_input/subtrace_daily_traded.RData')
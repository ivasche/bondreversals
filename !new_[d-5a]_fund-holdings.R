# ==================================================================================================
# [] FUND HOLDINGS PER BOND PER DAY (monthly, spread out)
#  - RUN TIME: ABOUT 20 MIN
# ==================================================================================================

library(data.table); library(xts); library(fastmatch); library(magrittr); library(pbapply)
source('./codes/delete.R')

# load bond data just to extract cusip numbers + bond-month obs
load('./data_input/subtrace.RData')
trace_cusips = unique(subtrace[, .(substr(cusip_id, 1, 8), cusip_id, trd_exctn_dt)])
names(trace_cusips)[1:2] = c('cusip8', 'cusip')
trace_cusips[, ym := as.yearmon(trd_exctn_dt)]
remove(submergent, subtrace)
gc()

# load funds data
load('./data_input/mf_holdings.RData')

# leave only what is needed
data = delete(data, keep.idxs = which(data$cusip %fin% unique(trace_cusips$cusip8)))
data[, hold := (nbr_shares > 0 | market_val > 0 | percent_tna > 0)]
data = delete(data, keep.idxs = which(data$hold == T))
data[, hold := as.numeric(hold)]
data[, c('percent_tna', 'nbr_shares', 'market_val', 'security_name', 'coupon') := NULL]
gc()

data[, report_dt := as.IDate(as.character(report_dt), format = '%Y%m%d')]
data[, ym := as.yearmon(report_dt)]
data[, report_dt := NULL]
remove(cusips, submergent)
gc()

## define a function to create a resulting dataset per each fund separately
gethold = function(x) {
  
  df = x[, seq(min(ym), max(ym), 1/12), cusip]
  names(df)[2] = 'ym'
  df[, hold := 1]
  
  return(df)
  
}

## apply
data <- split(data, f = data$crsp_portno)
data <- pblapply(data, gethold)
data <- rbindlist(data, use.names = TRUE, fill = TRUE, idcol = 'crsp_portno')

## extract no. funds per cusip per month
data <- data[, .N, .(cusip, ym)]
names(data)[3] = 'nfunds'

## add 9-digit cusip
names(trace_cusips)[1:2] = c('cusip', 'cusip_id')
trace_cusips <- data[trace_cusips, on = c('cusip', 'ym')]
trace_cusips[is.na(nfunds), nfunds := 0]
trace_cusips[, c('cusip', 'ym') := NULL]
gc()

trace_cusips <- trace_cusips[, .(cusip_id, trd_exctn_dt, nfunds)]
nfunds <- trace_cusips
remove(trace_cusips)
gc()

save(nfunds, file = './data_input/nfunds.RData')
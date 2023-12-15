# ==================================================================================================
# [1b] FILL OUTSTANDING AMOUNTS AND RATINGS (DAILY)
#  - RUN TIME: around 30 mins
# ==================================================================================================

# load packages
library(data.table); library(magrittr); library(fastmatch); library(pbapply); library(xts);
library(doParallel); library(iterators); library(timeDate); library(bizdays); library(fasttime)
library(xtable); library(stargazer); library(haven)

# set locale
Sys.setlocale(locale = 'US')

# load user-defined function
source('./codes/delete.R')
source('./codes/ratfunc.R')
source('./codes/desc.R')

# --------------------------------------------------------------------------------------------------
# read trace records, shrink to daily
# --------------------------------------------------------------------------------------------------

load('./data_input/daily_trace+ytm.RData')
subtrace = dtrace
remove(dtrace)
gc()

# --------------------------------------------------------------------------------------------------
# add and clean historical outstanding amounts
# --------------------------------------------------------------------------------------------------

## read changes in historical outstanding amounts
histao = read_sas('./data_input/fisd_amt_out_hist.sas7bdat')
histao = as.data.table(histao)
histao = histao[ISSUE_ID %in% submergent[, ISSUE_ID]]
histao = submergent[,.(ISSUE_ID, cusip_id)][histao, ,on = 'ISSUE_ID']

## pick what to add from this file
toadd = histao[, .(cusip_id, ACTION_TYPE, EFFECTIVE_DATE, AMOUNT_OUTSTANDING/1000)]
names(toadd) = c('cusip_id', 'action_type_hist', 'action_date_hist', 'size_hist')
toadd[, trd_exctn_dt := action_date_hist]

## also add amount at issuance (not always in histao!)
toadd2 = submergent[, .(cusip_id, ACTION_TYPE, OFFERING_DATE, OFFERING_AMT/1000)]
names(toadd2) = c('cusip_id', 'action_type_hist', 'action_date_hist', 'size_hist')
toadd2[, action_date_hist := as.Date(as.character(action_date_hist), format = '%Y%m%d')]
toadd2[, action_type_hist := 'I2']
toadd2[, trd_exctn_dt := action_date_hist]

## merge 2
toadd = rbind(toadd, toadd2)
remove(toadd2)

## clean giving priority to issuance from mergent master file
toadd[, mult := .N, .(cusip_id, trd_exctn_dt)]

tokeep = toadd[mult == 1]
toadd = toadd[mult > 1]

pickentry = function(x) {
  size = x$size_hist
  type = x$action_type_hist
  
  A = (type == 'I2')
  if (sum(A) == 0 & sd(A) == 0) {A[1] = T}
  if (sum(A) == 0 & sd(A) >0) {A[which.max(A)] = T}
  
  return(list(A))
}

toadd[, incl := pickentry(.SD), .(cusip_id, trd_exctn_dt)]
toadd = toadd[incl == T]
toadd[, incl := NULL]
toadd = rbind(tokeep, toadd)
toadd[, mult := NULL]
remove(histao, tokeep)

setkey(subtrace, cusip_id, trd_exctn_dt)
setkey(toadd, cusip_id, trd_exctn_dt)
gc()

subtrace = toadd[subtrace,, on = c('cusip_id', 'trd_exctn_dt'), roll = +Inf]
remove(toadd)
gc()

# about 1k of 18 mln bond-day obs have missing historical size, we treat this step by step
subtrace[, .N, is.na(size_hist)]

# flag bonds with no historical AO at all (2 issues)
nohistao = subtrace[, .(sum(!is.na(size_hist)), sum(is.na(size_hist))), cusip_id][V1 == 0 & V2 > 0][, cusip_id]
subtrace = delete(subtrace, which(!(subtrace$cusip_id %fin% nohistao)))
remove(nohistao)
gc()

# for the rest, fill size with issue size
fillna = function(x) {A = x; A[is.na(x)] = A[!is.na(x)][1]; return(A)}
subtrace[, size_hist := fillna(size_hist), .(cusip_id)]

# remove non-needed columns
subtrace[, c('action_type_hist', 'action_date_hist') := NULL]

# cut submergent to the subset of bonds
submergent = submergent[cusip_id %fin% unique(subtrace[, cusip_id])]

# --------------------------------------------------------------------------------------------------
# add historical ratings
# --------------------------------------------------------------------------------------------------

# read historical ratings, leave only issues that are in TRACE
histrat = fread('./data_input/ratings.csv')
histrat = histrat[ISSUE_ID %in% submergent[, ISSUE_ID]]

# add historical ratings
histrat[, RATING_DATE := as.Date(as.character(RATING_DATE), format = '%Y%m%d')]
histrat[, irow := .I]
histrat[, uni := rat(.SD), by = irow]
histrat[, uninum := ratnum(.SD), by = irow]

toadd = histrat[, .(COMPLETE_CUSIP, ISSUE_ID, RATING_DATE, uni, uninum)]
names(toadd) = c('cusip_id', 'issue_id', 'trd_exctn_dt', 'rating', 'rating_num')
remove(histrat)
gc()

# break ties in ratings
fillmultrat = function(x) { 
  
  if (length(x)  == 1) {
    
    return(1)
    
  } else if (sum(x == 23) == length(x)) {
    
    return(c(1, rep(0, length(x)-1)))
    
  } else {
    
    return(x == max(x[x<=22]))
    
  }

}

toadd[, tokeep := fillmultrat(rating_num), .(cusip_id, trd_exctn_dt)]
toadd = toadd[tokeep == 1]
toadd[, tokeep := NULL]
toadd = unique(toadd)

# sort before merging (critical!)
setkey(subtrace, cusip_id, trd_exctn_dt)
setkey(toadd, cusip_id, trd_exctn_dt)

# replcace numerical codes for missing ratings with NAs
toadd[rating == 'NR', rating := NA]
toadd[rating_num == 23, rating_num := NA]

# merge
subtrace = toadd[subtrace,, on = c('cusip_id', 'trd_exctn_dt'), roll = +Inf]
remove(toadd)
gc()

# read historical sp ratings
histrat = fread('./data_input/[data]_spratings_hist.csv')
histrat = histrat[, .(CUSIP9, ratingsymbol, ratingDate, IsCurrent)]
names(histrat) = c('cusip_id', 'RATING', 'trd_exctn_dt', 'iscurrent')
histrat[, trd_exctn_dt := as.Date(as.character(trd_exctn_dt), format = '%Y%m%d')]
histrat[, irow := .I]
histrat[, uni := rat(.SD), by = irow]
histrat[, uninum := ratnum(.SD), by = irow]
histrat = histrat[!cusip_id == '']
histrat[, c('RATING', 'iscurrent', 'irow') := NULL]
gc()

# break ties
histrat[, tokeep := fillmultrat(uninum), .(cusip_id, trd_exctn_dt)]
histrat = histrat[tokeep == 1]
histrat[, tokeep := NULL]
histrat = unique(histrat)

setkey(subtrace, cusip_id, trd_exctn_dt)
setkey(histrat, cusip_id, trd_exctn_dt)

# replcace numerical codes for missing ratings with NAs
histrat[uni == 'NR', uni := NA]
histrat[uninum == 23, uninum := NA]

# merge and add sp ratings where mergent ratings are absent
subtrace = histrat[subtrace,, on = c('cusip_id', 'trd_exctn_dt'), roll = +Inf]
subtrace[is.na(rating) & !is.na(uni), rating := uni]
subtrace[is.na(rating_num) & !is.na(uninum), rating_num := uninum]

subtrace[, c('uni', 'uninum') := NULL]
remove(histrat)
gc()

# read historical sp ratings
histrat_e = fread('./data_input/[data]_spentityratings_hist.csv')
histrat_e = histrat_e[, .(CUSIP6, ratingsymbol, ratingdate)]
names(histrat_e) = c('cusip6', 'RATING', 'trd_exctn_dt')
histrat_e[, trd_exctn_dt := as.Date(as.character(trd_exctn_dt), format = '%Y%m%d')]
histrat_e[, irow := .I]
histrat_e[, uni := rat(.SD), by = irow]
histrat_e[, uninum := ratnum(.SD), by = irow]
histrat_e[, c('RATING', 'irow') := NULL]
gc()

# break ties
histrat_e[, tokeep := fillmultrat(uninum), .(cusip6, trd_exctn_dt)]
histrat_e = histrat_e[tokeep == 1]
histrat_e[, tokeep := NULL]
histrat_e = unique(histrat_e)

# replcace numerical codes for missing ratings with NAs
histrat_e[uni == 'NR', uni := NA]
histrat_e[uninum == 23, uninum := NA]

# remove bonds with missing issuer
histrat_e = histrat_e[!cusip6==""]

# merge
subtrace[, cusip6 := substr(cusip_id, 1, 6)]

setkey(subtrace, cusip6, cusip_id, trd_exctn_dt)
setkey(histrat_e, cusip6, trd_exctn_dt)

subtrace = histrat_e[subtrace,, on = c('cusip6', 'trd_exctn_dt'), roll = +Inf]

subtrace[is.na(rating) & !is.na(uni), rating := uni]
subtrace[is.na(rating_num) & !is.na(uninum), rating_num := uninum]

subtrace[, c('uni', 'uninum') := NULL]
remove(histrat_e)
gc()

save(file = './data_input/subtrace_daily_outstanding_rating.RData', subtrace, submergent)
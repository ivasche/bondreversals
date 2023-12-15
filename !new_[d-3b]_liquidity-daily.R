# ==================================================================================================
# [3b] COMPUTE LIQUIDITY INDICATORS FOR DAILY DATA 
#  - RUN TIME: 2 HOURS
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
remove(subtrace, submergent)
gc()

# remove what is not needed
trace[, c('company_symbol') := NULL]
setkey(trace, cusip_id, trd_exctn_dt)
gc()

## mark small trades (second definition)
trace[, small100 := F]
trace[entrd_vol_qt <= 10000, small100 := T]

# [1] compute average bid-ask -----------------------------------------
A = trace[rpt_side_cd == 'B', mean(rptd_pr), .(cusip_id, trd_exctn_dt)]
names(A)[3] = 'bp'
B = trace[rpt_side_cd == 'S', mean(rptd_pr), .(cusip_id, trd_exctn_dt)]
names(B)[3] = 'sp'
avg_ba = merge(A, B, all = T)
avg_ba[, avg_ba := 100*(sp-bp)/(0.5*sp+0.5*bp)]
avg_ba[, bp := NULL]; avg_ba[, sp := NULL]
remove(A, B)

# clean average bid-ask
avg_ba[avg_ba < 0, avg_ba := NA]
avg_ba[, avg_ba := winsorize(avg_ba, 0, 0.975)]

# [2] compute average bid-ask for large trades ------------------------
A = trace[rpt_side_cd == 'B' & small100 == F, mean(rptd_pr), .(cusip_id, trd_exctn_dt)]
names(A)[3] = 'bp'
B = trace[rpt_side_cd == 'S' & small100 == F, mean(rptd_pr), .(cusip_id, trd_exctn_dt)]
names(B)[3] = 'sp'
avg_ba_large = merge(A, B, all = T)
avg_ba_large[, avg_ba_lrg := 100*(sp-bp)/(0.5*sp+0.5*bp)]
avg_ba_large[, bp := NULL]
avg_ba_large[, sp := NULL]
remove(A, B)

# clean average bid-ask 
avg_ba_large[avg_ba_lrg < 0, avg_ba_lrg := NA]
avg_ba_large[, avg_ba_lrg := winsorize(avg_ba_lrg, 0, 0.975)]
gc()

# [3] compute gamma --------------------------------------------------
gamma = trace[, -var(log(rptd_pr)-shift(log(rptd_pr)),
                     shift(log(rptd_pr))-shift(log(rptd_pr),2),
                     na.rm = T), .(cusip_id, trd_exctn_dt)]
names(gamma)[3] = 'gamma'
N = trace[, .(.N), .(cusip_id, trd_exctn_dt)]
gamma = gamma[N,,on = c('cusip_id', 'trd_exctn_dt')]
remove(N)
gamma[, gamma := gamma*10^4]

# clean gamma
gamma[N<=4, gamma := NA]
gamma[, gamma := winsorize(gamma, 0.001, 0.975)]
gamma[, N := NULL]

# [4] compute roll ---------------------------------------------------
roll = trace[, -var(100*(rptd_pr/shift(rptd_pr)-1),
                    shift(100*(rptd_pr/shift(rptd_pr)-1)),
                    na.rm = T), .(cusip_id, trd_exctn_dt)]
names(roll)[3] = 'cov'
roll[, roll := 2*sqrt(cov*as.numeric(cov>0))]
roll[, cov := NULL]
N = trace[, .(.N), .(cusip_id, trd_exctn_dt)]
roll = roll[N,,on = c('cusip_id', 'trd_exctn_dt')]
remove(N)

# clean roll
roll[N<=4, roll := NA]
roll[, roll := winsorize(roll, 0.001, 0.975)]
roll[, N := NULL]
gc()

# [5] compute iqr --------------------------------------------------
iqr = trace[, 100*(quantile(rptd_pr, 0.75) - quantile(rptd_pr, 0.25))/mean(rptd_pr),
            .(cusip_id, trd_exctn_dt)]
names(iqr)[3] = 'iqr'
N = trace[, .(.N), .(cusip_id, trd_exctn_dt)]
iqr = iqr[N,,on = c('cusip_id', 'trd_exctn_dt')]
remove(N)

# clean iqr
iqr[N<=2, iqr := NA]
iqr[, iqr := winsorize(iqr, 0, 0.975)]
iqr[, N := NULL]

# remove what is not needed
remove(trace)
gc()

# merge all
dliq = merge(merge(merge(avg_ba, iqr, all = T), merge(gamma, roll, all = T), all = T), avg_ba_large, all = T)
remove(avg_ba_large, avg_ba, roll, iqr, gamma)
gc()

dliq[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]

# save
save(dliq, file = './data_input/liquidity_daily.RData')
# =============================================================================================
# - THIS CODE COMPUTES DAILY ACCRUED INTEREST AND YIELDS FOR ALL SAMPLE BONDS OBSERVED DAILY
# - DAILY PRICES ARE VWAPs. ABOUT 9 HOURS
# =============================================================================================

# =============================================================================================
# (0) Load prerequisites and read clean tick data
# =============================================================================================

# load packages
library(data.table); library(magrittr); library(fastmatch); library(pbapply); library(xts)
library(doParallel); library(iterators); library(timeDate); library(bizdays); library(fasttime)
library(xtable); library(stargazer); library(haven); library(nleqslv); library(BondValuation)
library(bizdays); library(Rcpp); library(RcppArmadillo)

# set locale
Sys.setlocale(locale = 'US')

# load user-defined function
source('./codes/delete.R')
source('./codes/ratfunc.R')
source('./codes/desc.R')
sourceCpp('./codes/ekfrun_cpp_singleobs.cpp')
source('./codes/[func]_CIR2+D99_KFs.R')

# define two thresholds: how many trading days per bond we need to observe (numdays_curtoff)
# and how many overlapping trading days for a pair of callable/non-callable bonds  
# of the same issuer we need to observe (numobs_cutoff)

#numdays_cutoff = 36
#numobs_cutoff  = 24

# read tick data
load('./data_input/subtrace_tick.RData')

# extract bond-day observations and remove the rest for now
subtrace[, `:=` (vwap = weighted.mean(x = rptd_pr, w = entrd_vol_qt)), .(cusip_id, trd_exctn_dt)]
dtrace = unique(subtrace[, .(cusip_id, trd_exctn_dt, vwap)])
remove(subtrace)
gc()

# leave only bonds with at least numdays_cutoff daily observations
#tokeep = dtrace[, .N, cusip_id][N >= numdays_cutoff][, cusip_id]
#dtrace = dtrace[cusip_id %fin% tokeep]
#remove(tokeep)
#gc()

# cut submergent further
submergent = submergent[cusip_id %fin% unique(dtrace[, cusip_id])]

# =============================================================================================
# (2) Add to TRACE bond-day panel whatever is useful from MERGENT datasets
# =============================================================================================

# add what we know about callability from mergent masterfile
# scaling of 'amounts' means to bring them all to mln USD throught this file
toadd = submergent[, .(cusip_id, ISSUE_ID, ISSUER_ID, REDEEMABLE, OFFERING_DATE, DATED_DATE, FIRST_INTEREST_DATE,
                       LAST_INTEREST_DATE, MATURITY, OFFERING_AMT/1000, COUPON,
                       INTEREST_FREQUENCY)]
names(toadd) = c('cusip_id', 'ISSUE_ID', 'issuer', 'callable', 'issue_date', 'dated_date', 'first_int_date',
                 'last_int_date', 'maturity', 'issue_size', 'coupon', 'int_freq')
dtrace = toadd[dtrace,, on = 'cusip_id']
remove(toadd)
gc()

# change date fields to date format
dtrace[, issue_date := as.IDate(as.character(issue_date), format = '%Y%m%d')]
dtrace[, dated_date := as.IDate(as.character(dated_date), format = '%Y%m%d')]
dtrace[, first_int_date := as.IDate(as.character(first_int_date), format = '%Y%m%d')]
dtrace[, last_int_date := as.IDate(as.character(last_int_date), format = '%Y%m%d')]
dtrace[, maturity := as.IDate(as.character(maturity), format = '%Y%m%d')]
dtrace[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]
gc()

# remove trades on maturity date
dtrace = dtrace[!trd_exctn_dt == maturity]
submergent = submergent[cusip_id %fin% unique(dtrace[, cusip_id])]
gc()

# ==========================================================================================================================
# (2) Calculate accrued interest
# ==========================================================================================================================

# pick bond characteristics
toadd = submergent[, .(cusip_id, MATURITY, OFFERING_DATE, DATED_DATE, FIRST_INTEREST_DATE, LAST_INTEREST_DATE, COUPON,
                       INTEREST_FREQUENCY)]
datecols = c("MATURITY", "OFFERING_DATE", "DATED_DATE", "FIRST_INTEREST_DATE", "LAST_INTEREST_DATE")
toadd[, (datecols) := lapply(.SD, function(x) {as.IDate(as.character(x), format = '%Y%m%d')}),, .SDcols = datecols]

# replace missing offering dates with dated dates (usual assumption)
toadd[is.na(OFFERING_DATE), OFFERING_DATE := DATED_DATE]

# pick only non-zero-coupon bonds
toadd = toadd[!(COUPON == 0)]

# index rows
toadd[, irow := .I]

# bonds with missing maturity but non-missing last interest date (likely error reporting)
toadd[is.na(MATURITY) & !is.na(LAST_INTEREST_DATE), MATURITY := seq(from = LAST_INTEREST_DATE, length.out = 2,
                                                                    by = paste(12/INTEREST_FREQUENCY, 'month'))[2], irow]

# if first interest date == maturity, then it also equals the last interest date
toadd[!is.na(MATURITY) & is.na(LAST_INTEREST_DATE) & (FIRST_INTEREST_DATE == MATURITY),
      LAST_INTEREST_DATE := FIRST_INTEREST_DATE]

# if last interest date still missing (with non-na maturity) then fill with the seq from first to maturity
toadd[!is.na(MATURITY) & is.na(LAST_INTEREST_DATE),
      LAST_INTEREST_DATE := last(seq(from = FIRST_INTEREST_DATE, to = MATURITY, by = paste(12/INTEREST_FREQUENCY, 'month'))),
      irow]

# pure perpetuity
prptl = toadd[is.na(MATURITY), cusip_id]
toadd[cusip_id %in% prptl, MATURITY := seq(from = OFFERING_DATE, length.out = 2, by = '100 year')[2], irow]
toadd[cusip_id %in% prptl, LAST_INTEREST_DATE := seq(from = MATURITY, length.out = 2,
                                                     by = paste(12/INTEREST_FREQUENCY, 'month'))[2], irow]

# DO SMTH with ODD coupons
toadd[, ODD := last(seq(from = FIRST_INTEREST_DATE, to = LAST_INTEREST_DATE, by = paste(12/INTEREST_FREQUENCY, 'month'))),
      irow]
toadd[, ODD := !(as.numeric(ODD - LAST_INTEREST_DATE) == 0)]

# create coupon schedules FOR COUPON BONDS ONLY
L = pblapply(as.list(1:nrow(toadd)), function(i) {
  
  A = AnnivDates(Em = toadd$OFFERING_DATE[i], Mat = toadd$MATURITY[i], CpY = toadd$INTEREST_FREQUENCY[i],
                 FIPD = toadd$FIRST_INTEREST_DATE[i], LIPD = toadd$LAST_INTEREST_DATE[i],
                 FIAD = toadd$DATED_DATE[i], RV = 100, Coup = toadd$COUPON[i], DCC = 9, EOM = F)$PaySched
  A[nrow(A), 2] = A[nrow(A), 2] + 100
  
  B = data.table(cusip_id = toadd$cusip_id[i],
                 dates = A[,1],
                 coupon = A[,2])
  return(B)
  
})

L = rbindlist(L)
L[, trd_exctn_dt := dates]

## add coupons (yearly) and no. payments per year + dated and maturity dates
L = toadd[, .(cusip_id, COUPON, INTEREST_FREQUENCY, OFFERING_DATE, DATED_DATE, FIRST_INTEREST_DATE, LAST_INTEREST_DATE,
              MATURITY, ODD)][L,,on = c('cusip_id')]
names(L)[2:9] = c('c_rate', 'int_freq', 'issue_date', 'dated_date', 'first_date', 'last_date', 'mat', 'odd')

# select unique bond-trading day pairs (only COUPON BONDS)
aint = dtrace[, .(cusip_id, trd_exctn_dt, vwap)]
aint[, ind := cusip_id]

# modify an accued interest function to return a plain numeric and vectorize this
numAccrInt = function(StartDate, EndDate, Coup, DCC, CpY, RV, YearNCP, EOM) {
  AccrInt(StartDate=StartDate, EndDate=EndDate, Coup=Coup, DCC=DCC, CpY=CpY, RV=RV, YearNCP=YearNCP, EOM=EOM)[[1]]
}
v_AccrInt = Vectorize(numAccrInt,vectorize.args = c('StartDate', 'EndDate', 'Coup', 'CpY', 'YearNCP'))

## define a function that creates a list (as many elements as there're trading dates) with times to payments and amounts
cashflows = function(y, sched) {
  
  B = as.data.table(sched)
  B[, tomat := as.numeric(dates - y)/360]
  
  list(tomat = B[tomat>0, tomat], cfs = B[tomat>0, coupon])
  
}

## read risk-free curve
yrf = read.table('./data_input/YC.txt')
yrf = cbind(as.IDate(rownames(yrf)), yrf)
yrf = as.data.table(yrf)
names(yrf)[1] = 'trd_exctn_dt'

## make a sequence of dates (incl weekends)
dts = seq.Date(from = yrf$trd_exctn_dt[1], to = yrf$trd_exctn_dt[nrow(yrf)], by = 'days')
dts = as.data.table(dts)
names(dts) = 'trd_exctn_dt'

## carry last obs forward (fill weekends)
yrf = yrf[dts,, on = 'trd_exctn_dt', roll = +Inf]
remove(dts)
gc()

## extend to the max maturity
maxmat_trace = 1300
maxmat_yrf = as.numeric(substr(names(yrf)[ncol(yrf)], 2, 4))
cnms = paste0("M", seq(maxmat_yrf + 1, maxmat_trace, 1))

yrf = cbind(yrf, as.data.table(matrix(rep(yrf[, ncol(yrf)], length(cnms)), ncol = length(cnms))))
names(yrf)[which(grepl('^V', names(yrf)))] = cnms

## define a funciton that gives a list of yields, accrued interests and coupons per bond for all trading days
coupfn = function(x) {
  
  # compute accrued interest and coupons
  schedule = L[cusip_id == x$cusip_id[1]]
  nrowsch = nrow(schedule)
  
  D1 = schedule[x, .(cusip_id, trd_exctn_dt, dates), on = c('cusip_id','trd_exctn_dt'),
                roll = +Inf]
  names(D1)[3] = 'pdate'
  
  D2 = schedule[x, .(cusip_id, trd_exctn_dt, dates, coupon, c_rate, int_freq, issue_date, dated_date, first_date,
                     last_date, mat),
                on = c('cusip_id','trd_exctn_dt'),
                roll = -Inf]
  names(D2)[3] = 'ndate'
  
  D = D1[D2,,on = c('cusip_id','trd_exctn_dt')]
  D = D[x[, .(cusip_id, trd_exctn_dt, vwap)],, on = c('cusip_id', 'trd_exctn_dt')]
  
  D[trd_exctn_dt <= first_date, pdate := dated_date]
  D[is.na(pdate) & trd_exctn_dt > first_date & trd_exctn_dt < ndate, pdate := first_date]
  
  D[, ai := v_AccrInt(StartDate = pdate, EndDate = trd_exctn_dt, Coup = c_rate, DCC = 9, CpY = int_freq, RV = 100,
                      YearNCP = year(ndate), EOM = F)]
  
  # create payment schedules for this bond
  bond.L = lapply(as.list(x$trd_exctn_dt), function(z) {cashflows(y = z, sched = schedule)})
  names(bond.L) = as.yearmon(x$trd_exctn_dt)
  
  # calculate risk-free prices
  yrf_in = yrf[trd_exctn_dt %in% D$trd_exctn_dt]
  
  dfs = lapply(as.list(seq_along(bond.L)), function(i) {
    
    list(dfs = exp(-as.numeric(yrf_in[i, paste0('M', round(12*bond.L[[i]]$tomat, 0)),
                                      with = F])/100*bond.L[[i]]$tomat))
    
    
  })
  
  bond.L = mapply(FUN = c, bond.L, dfs, SIMPLIFY = F)
  prf = as.numeric(sapply(bond.L, function(j) { cb_price_dfs_cpp(coup = j$cfs, dfs = j$dfs, notional = 1) }))
  
  # calculate yields
  yld = lapply(as.list(seq_along(D$trd_exctn_dt)), function(t) {
    
    list(ytm = 100*ytm_subcpp(coup = bond.L[[t]]$cfs,
                              coup_dates = bond.L[[t]]$tomat,
                              coup_freq = schedule$int_freq[1],
                              price = D$vwap[t]+D$ai[t],
                              notional = 1),
         ytm_rf = 100*ytm_subcpp(coup = bond.L[[t]]$cfs,
                                 coup_dates = bond.L[[t]]$tomat,
                                 coup_freq = schedule$int_freq[1],
                                 price = prf[t],
                                 notional = 1))
         })
  yld = rbindlist(yld)
  D = cbind(D, yld)
  
  # add maturity
  D[, mat := as.yearmon(mat) - as.yearmon(trd_exctn_dt)]
  D[!(trd_exctn_dt == pdate & trd_exctn_dt == ndate), coupon := 0]
  
  return(D[, .(cusip_id, trd_exctn_dt, ai, ytm, ytm_rf, coupon, mat)])
  
}

# run this thing 
aint_coup = aint[cusip_id %fin% toadd$cusip_id]
aint_coup = split(aint_coup, by = 'cusip_id')
res = pblapply(aint_coup, coupfn)
res = rbindlist(res)
gc()

# stuck results with the original file
aint = res[aint,, on = c('cusip_id', 'trd_exctn_dt')]
aint[, ind := NULL]

# add info on zero-coupon bonds
zc_cusips = unique(aint[is.na(mat), cusip_id])

toadd_zc = submergent[cusip_id %fin% zc_cusips, .(cusip_id, MATURITY)]
toadd_zc[, MATURITY := as.IDate(as.character(MATURITY), format = '%Y%m%d')]
names(toadd_zc)[2] = c('mat_zc')

aint = toadd_zc[aint,,on = c('cusip_id')]
aint[(cusip_id %fin% zc_cusips) & (mat_zc == trd_exctn_dt), coupon := 100]
aint[(cusip_id %fin% zc_cusips) & !(mat_zc == trd_exctn_dt), coupon := 0]
aint[cusip_id %fin% zc_cusips, ai := 0]

aint[cusip_id %fin% zc_cusips, mat := as.yearmon(mat_zc) - as.yearmon(trd_exctn_dt)]
aint[, mat_zc := NULL]

aint[cusip_id %fin% zc_cusips, ytm := 2*((100/(vwap+ai))^(1/(2*mat))-1)]

# clean 
remove(aint_coup, res, toadd, toadd_zc, datecols, prptl)
gc()

## add ytm_rf of zc bonds (same as risk-free rate with maturity T)
zc_mats = aint[cusip_id %fin% zc_cusips, round(12*mat)]

yrf_sub = as.data.table(aint[cusip_id %fin% zc_cusips, trd_exctn_dt])
names(yrf_sub) = 'trd_exctn_dt'
yrf_sub = yrf[yrf_sub,,on = 'trd_exctn_dt']

zc_yrfs = sapply(as.list(seq_along(zc_mats)), function(i) {as.numeric(yrf_sub[i, zc_mats[i]+2, with = F])})
aint[cusip_id %fin% zc_cusips, ytm_rf := zc_yrfs]

remove(zc_mats, yrf_sub, zc_yrfs, zc_cusips, L, yrf, cnms)
gc()

# add accrued interest to subtrace
dtrace = aint[dtrace,, on = c('cusip_id', 'trd_exctn_dt')]

# leave in master file only securities from trace 
submergent = submergent[cusip_id %in% unique(dtrace[, cusip_id])]
gc()

# clean
dtrace[, c('i.coupon', 'i.vwap') := NULL]

# save
save(dtrace, submergent, file = './data_input/daily_trace+ytm.RData')
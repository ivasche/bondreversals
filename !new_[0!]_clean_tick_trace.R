# ==================================================================================================
# [0] CLEAN TRACE, SUBMERGENT, and COMPUTE ACCRUED INTEREST
#  - RUN TIME: around 1 hour
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch); library(jrvFinance)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(doParallel); library(timeDate); library(bizdays); library(fasttime); library(haven)
library(BondValuation)

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')

# --------------------------------------------------------------------------------------------------
# make subsample of senior, fixed coupon, non-convertible, non-asset-backed debt
# --------------------------------------------------------------------------------------------------

# read mergent master file to match characteristics of securities in TRACE
mergent  = fread('./data_input/mergent_fisd_bond_issues.csv')

# clean the master file
submergent = mergent[FOREIGN_CURRENCY == 'N'                                    &
#                       RULE_144A == 'N'                                         &
                       PRIVATE_PLACEMENT == 'N'                                 &
                       YANKEE == 'N'                                            &
                       CANADIAN == 'N'                                          &  
                       BOND_TYPE %in% c('CMTN', 'CDEB', 'CMTZ', 'CZ', 'RNT',
                                        'USBN', 'PS', 'UCID')                   &
                       CONVERTIBLE == 'N'                                       &
                       ASSET_BACKED == 'N'                                      &
                       COUPON_TYPE %in% c('F','Z')                              &
                       ENHANCEMENT == 'N'                                       &
                       !SECURITY_PLEDGE == 'M'                                  
]

# fill zero coupons for zero-frequency payments
submergent[COUPON_TYPE == 'Z', COUPON := 0L]
submergent[COUPON_TYPE == 'Z', INTEREST_FREQUENCY := 0L]
submergent[COUPON == 0 | COUPON == 0L, INTEREST_FREQUENCY := 0L]
submergent[INTEREST_FREQUENCY == 0 & COUPON == 0, COUPON_TYPE := 'Z']

# remove fixed-coupon bonds with NA coupons
submergent = submergent[!(COUPON_TYPE == 'F' & is.na(COUPON))]

# replace NA interest frequency with fixed non-zero coupon with 22
submergent[is.na(INTEREST_FREQUENCY) & COUPON > 0 & COUPON_TYPE == 'F', INTEREST_FREQUENCY := 2L]

# remove fixed-coupon bonds with non-zero coupon but zero interest frequency (about 200 names of 315k)
# might be bullet repayment but this is very exotic
submergent = submergent[!(INTEREST_FREQUENCY == 0L & COUPON>0 & COUPON_TYPE == 'F')]

# assume that all weird interest frequencies == 2
submergent[INTEREST_FREQUENCY %fin% c(-1, 99), INTEREST_FREQUENCY := 2]

# for bonds with missing maturity (perpetuities), set maturity to OFFERING_DATE + 100 years
submergent[, irow := .I]
submergent = submergent[is.na(MATURITY), MATURITY := seq(from = as.IDate(as.character(OFFERING_DATE), format = '%Y%m%d'),
                                                          by = '100 year', length.out = 2)[2],
                        irow]

# --------------------------------------------------------------------------------------------------
# read trace records
# --------------------------------------------------------------------------------------------------

# read chunks and cut TRACE sample to only securities we indentified in MERGENT
load('./data_input/trace_04_08_new.RData')
trace_04_08 <- trace_04_08[trace_04_08$cusip_id %in% submergent$COMPLETE_CUSIP,]
gc()

load('./data_input/trace_09_12_new.RData')
trace_09_12 <- trace_09_12[trace_09_12$cusip_id %in% submergent$COMPLETE_CUSIP,]
gc()

load('./data_input/trace_12_14_new.RData')
trace_12_14 <- trace_12_14[trace_12_14$cusip_id %in% submergent$COMPLETE_CUSIP,]
gc()

load('./data_input/trace_15_17_new.RData')
trace_15_17 <- trace_15_17[trace_15_17$cusip_id %in% submergent$COMPLETE_CUSIP,]
gc()

load('./data_input/trace_17_19_new.RData')
trace_17_19 <- trace_17_19[trace_17_19$cusip_id %in% submergent$COMPLETE_CUSIP,]
gc()

load('./data_input/trace_19_20_new.RData')
trace_19_20 <- trace_19_20[trace_19_20$cusip_id %in% submergent$COMPLETE_CUSIP,]
gc()

traceall <- list(trace_04_08, trace_09_12, trace_12_14, trace_15_17, trace_17_19, trace_19_20)
remove(trace_04_08, trace_09_12, trace_12_14, trace_15_17, trace_17_19, trace_19_20)
gc()

traceall[[1]][, c('bond_sym_id', 'trd_exctn_tm') := NULL]
traceall[[2]][, c('bond_sym_id', 'trd_exctn_tm') := NULL]

traceall[[3]][, c('bond_sym_id', 'trd_exctn_tm', 'days_to_sttl_ct', 'agu_qsr_id', 'lckd_in_ind') := NULL]
traceall[[4]][, c('bond_sym_id', 'trd_exctn_tm', 'days_to_sttl_ct', 'agu_qsr_id', 'lckd_in_ind') := NULL]
traceall[[5]][, c('bond_sym_id', 'trd_exctn_tm', 'days_to_sttl_ct', 'agu_qsr_id', 'lckd_in_ind') := NULL]
traceall[[6]][, c('bond_sym_id', 'trd_exctn_tm', 'days_to_sttl_ct', 'agu_qsr_id', 'lckd_in_ind') := NULL]
gc()

subtrace <- rbindlist(traceall)
remove(traceall)
gc()

# --------------------------------------------------------------------------------------------------
# clean trace records
# --------------------------------------------------------------------------------------------------

# leave in master file only securities from trace and fields defined above
submergent <- submergent[COMPLETE_CUSIP %in% subtrace[,.(.N),cusip_id]$cusip_id,]

# remove not needed stuff
remove(mergent)
gc()

# rename cusip_id part in master file
names(submergent)[grep('COMPLETE', names(submergent))] = 'cusip_id'

## remove cusips with a single trading day
subtrace <- delete(subtrace,
                  !(subtrace$cusip_id %in% subtrace[,length(unique(trd_exctn_dt)),
                                                    cusip_id][V1 == 1][,cusip_id]))
gc()

# remove bonds traded at <5 or >1000
t5 <- subtrace[, .(cusip_id, rptd_pr)][, .(sum(rptd_pr <= 5), sum(rptd_pr >= 1000), .N), cusip_id]
t5[, V1 := V1/N*100]; t5[, V2 := V2/N*100]

subtrace[, ind := F]
subtrace[cusip_id %in% t5[V1 == 100 | V2 == 100, cusip_id], ind := T]
subtrace <- delete(subtrace, which(subtrace$ind == FALSE))
remove(t5); subtrace[, ind := NULL]
gc()

# remove trades that are < 10k USD
#subtrace[, ind := F]
#subtrace[entrd_vol_qt < 10000, ind := T]
#subtrace = delete(subtrace, which(subtrace$ind == FALSE))
#subtrace[, ind := NULL]
#gc()

# remove not needed stuff
subtrace[, c('rptg_party_type') := NULL]
gc()

# add dated_date and maturity
toadd <- submergent[, .(cusip_id, MATURITY, DATED_DATE)]
subtrace <- toadd[subtrace, , on='cusip_id']
remove(toadd)
gc()

# remove trades before dated_date and after maturity
subtrace <- subtrace[trd_exctn_dt<=MATURITY & trd_exctn_dt>=DATED_DATE]
subtrace[, c('MATURITY', 'DATED_DATE') := NULL]
gc()

# save submergent separately
#save(submergent, file = './data_input/submergent_tick.RData')
#remove(submergent)
#gc()

# convert date field to date
#subtrace[, trd_exctn_dt := as.IDate(as.character(trd_exctn_dt), format = '%Y%m%d')]
#gc()

## save
save(subtrace, submergent, file = './data_input/subtrace_tick.RData')

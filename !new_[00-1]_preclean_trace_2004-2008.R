# ==================================================================================================
# [0-1] THIS PIECE CLEANS 2004-2008 DATA
#  - ABOUT 15 MINUTES  
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch); library(utils)

# set locale
Sys.setlocale(locale = 'US')

# load delete-by-reference code
source('./codes/delete.R')

# --------------------------------------------------------------------------------------------------
# Do Dick-Nielsen filtering for 2004-2008
# --------------------------------------------------------------------------------------------------

## read data
trace12 = fread(input = './data_input/trace_2004_2008.csv', header = T, stringsAsFactors = F)
trace12[,irow := .I]

## delete columns with unapplied filters
trace12[, c('stlmnt_dt', 'yld_pt', 'yld_sign_cd', 'sub_prdct', 'trd_mod_3', 'trd_mod_4',
            'sale_cndtn2_cd', 'lckd_in_ind', 'dissem_fl', 'pr_trd_dt', 'bloomberg_identifier',
            'ats_indicator', 'first_trade_ctrl_num', 'first_trade_ctrl_date') := NULL]

## index cancellations + original trades, corrections
trace12_out1 = trace12[trc_st %in% c('C','Y')]
out1row = c(trace12_out1$irow,
            trace12[trace12_out1, irow, on = c('trd_rpt_dt', msg_seq_nb = 'orig_msg_seq_nb')])
out2row = trace12[trc_st == 'C', irow]
remove(trace12_out1)

## index all trades to keep
ind = !fmatch(trace12$irow, c(out1row, out2row), nomatch = 0)>0
remove(out1row, out2row)

## remove the rest
trace12 = delete(trace12, ind)
remove(ind)

## index reversals
trace12_out2 = trace12[asof_cd == 'R']
ind = which(!trace12$asof_cd == 'R')
trace12      = delete(trace12, ind)

outrow = trace12[trace12_out2,
                 on = c('trd_exctn_dt', 'cusip_id', 'trd_exctn_tm', 'rptd_pr', 'entrd_vol_qt',
                        'rpt_side_cd', 'cntra_mp_id')]
outrow = outrow[trd_exctn_dt<i.trd_rpt_dt, .(irow, i.irow)]
outrow = unique(outrow, by = 'i.irow')

## index all trades to keep
ind = !fmatch(trace12$irow, outrow$irow, nomatch = 0)>0

## remove the rest
trace12 = delete(trace12, ind)
trace12[, irow := NULL]
remove(ind, outrow, trace12_out2)

## delete message sequence number column and reporting times
trace12[, c('orig_msg_seq_nb', 'trd_rpt_dt', 'trd_rpt_tm') := NULL]
sapply(trace12, function(x) {sum(is.na(x))})

# --------------------------------------------------------------------------------------------------
# Remove agency transactions without comissions
# --------------------------------------------------------------------------------------------------
trace12[, agency := buy_cpcty_cd]
trace12[rpt_side_cd == 'S', agency := sell_cpcty_cd]
trace12[, irow := .I]
ind = trace12[!(agency == 'A' & cntra_mp_id == 'C' & cmsn_trd == 'N'), irow]
trace12 = delete(trace12, ind)
trace12[, irow := NULL]
trace12[, agency := NULL]
#trace12[,c('cmsn_trd') := NULL]

# --------------------------------------------------------------------------------------------------
# apply filters
# --------------------------------------------------------------------------------------------------

# regular filters
trace12[,irow := .I]

trace12[cusip_id == '', irow := NA]
trace12[is.na(rptd_pr), irow := NA]

trace12[wis_fl == 'Y', irow := NA]
trace12[trdg_mkt_cd %in% c('S2', 'P1', 'P2'), irow := NA]
trace12[spcl_trd_fl == 'Y', irow := NA]
trace12[scrty_type_cd == 'E', irow := NA]
trace12[cmsn_trd == 'Y', irow := NA]
trace12[days_to_sttl_ct > 6, irow := NA]
trace12[agu_qsr_id == 'G', irow := NA]
trace12[sale_cndtn_cd == 'C', irow := NA]
trace12[trc_st %in% c('X', 'R'), irow := NA]
trace12 = delete(trace12, keep.idxs = trace12[!is.na(irow), irow])

# inter-dealer trades
trace12[,irow := .I]
trace12[cntra_mp_id == 'D' & rpt_side_cd == 'B', irow := NA]
trace12[cntra_mp_id == 'D' & rpt_side_cd == 'S', rpt_side_cd := 'D']
trace12 = delete(trace12, keep.idxs = trace12[!is.na(irow), irow])

trace12[, c('irow', 'wis_fl', 'trdg_mkt_cd', 'spcl_trd_fl', 'scrty_type_cd', 'sale_cndtn_cd',
            'agu_qsr_id', 'days_to_sttl_ct', 'msg_seq_nb', 'asof_cd', 'trc_st',
            'cmsn_trd', 'buy_cmsn_rt', 'buy_cpcty_cd', 'sell_cmsn_rt', 'sell_cpcty_cd',
            'cntra_mp_id') := NULL]

## remove some commas from bond names
trace12[grep(',',trace12$bond_sym_id), bond_sym_id := sub(',', '.', bond_sym_id)]
trace12[grep(',',trace12$company_symbol), company_symbol := sub(',', '.', company_symbol)]

## save file
#fwrite(trace12, './data_input/trace_2004-2008_clean.csv')
trace_04_08 = trace12
remove(trace12)
save(trace_04_08, file= './data_input/trace_04_08.RData')
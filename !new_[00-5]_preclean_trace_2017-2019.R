# ==================================================================================================
# [0-5] THIS PIECE CLEANS 2017-2019 DATA
#  - ABOUT 15 MINUTES  
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)

# set locale
Sys.setlocale(locale = 'US')

# load delete-by-reference code
source('./codes/delete.R')

# --------------------------------------------------------------------------------------------------
# Do Dick-Nielsen filtering for 2017-2019
# --------------------------------------------------------------------------------------------------

## read data
trace12 = fread(input = './data_input/trace_2017-2019.csv', header = T, stringsAsFactors = F)
trace12[,irow := .I]

## delete columns with unapplied filters
trace12[, c('stlmnt_dt', 'yld_pt', 'yld_sign_cd', 'sub_prdct', 'trd_mod_3', 'trd_mod_4',
            'sale_cndtn2_cd', 'lckd_in_ind', 'dissem_fl', 'pr_trd_dt', 'bloomberg_identifier',
            'ats_indicator', 'first_trade_ctrl_num', 'first_trade_ctrl_date') := NULL]

## index cancellations and corrections + original trades
trace12_out1 = trace12[trc_st %in% c('X','C')]
out1row = trace12[trace12_out1, irow, on = c('cusip_id', 'entrd_vol_qt', 'rptd_pr', 'trd_exctn_dt',
                                             'trd_exctn_tm', 'rpt_side_cd', 'cntra_mp_id',
                                             'msg_seq_nb')]

## index reversals + original trades
trace12_out2 = trace12[trc_st %in% c('Y')]
out2row = trace12[trace12_out2, irow, on = c('cusip_id', 'entrd_vol_qt', 'rptd_pr', 'trd_exctn_dt',
                                             'trd_exctn_tm', 'rpt_side_cd', 'cntra_mp_id',
                                             msg_seq_nb = 'orig_msg_seq_nb')]
out2row = c(out2row, trace12_out2[,irow])

## index all trades to keep
ind = !fmatch(trace12$irow, c(out1row, out2row), nomatch = 0)>0
remove(trace12_out2, trace12_out1, out2row, out1row)

## remove the rest
trace12 = delete(trace12, ind)
trace12[, irow := NULL]
remove(ind)

## delete message sequence number column and reporting times
trace12[, c('orig_msg_seq_nb', 'trd_rpt_dt', 'trd_rpt_tm', 'agu_qsr_id') := NULL]
trace12[, c('scrty_type_cd', 'sale_cndtn_cd') := NULL]
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
trace12[cmsn_trd == 'Y', irow := NA]
trace12[trc_st %in% c('X', 'R'), irow := NA]

trace12 = delete(trace12, keep.idxs = trace12[!is.na(irow), irow])
gc()

# inter-dealer trades
trace12[,irow := .I]
trace12[cntra_mp_id == 'D' & rpt_side_cd == 'B', irow := NA]
trace12[cntra_mp_id == 'D' & rpt_side_cd == 'S', rpt_side_cd := 'D']
trace12 = delete(trace12, keep.idxs = trace12[!is.na(irow), irow])

trace12[, c('irow', 'wis_fl', 'trdg_mkt_cd', 'spcl_trd_fl', 'buy_cmsn_rt', 'sell_cmsn_rt',
            'cmsn_trd', 'msg_seq_nb', 'asof_cd', 'trc_st', 'buy_cpcty_cd', 'sell_cpcty_cd',
            'cntra_mp_id') := NULL]

trace12[grep(',',trace12$bond_sym_id), bond_sym_id := sub(',', '.', bond_sym_id)]
trace12[grep(',',trace12$company_symbol), company_symbol := sub(',', '.', company_symbol)]

## save file
#fwrite(trace12, './data_input/trace_2017-2019_clean.csv')
trace_17_19 = trace12
remove(trace12)
save(trace_17_19, file= './data_input/trace_17_19.RData')
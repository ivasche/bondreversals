delete <- function(DT, keep.idxs){
  cols <- copy(names(DT))
  DT_subset <- DT[[1]][keep.idxs] %>% as.data.table
  setnames(DT_subset, ".", cols[1])
  for (col in cols){
    DT_subset[, (col) := DT[[col]][keep.idxs]]
    set(DT, NULL, col, NULL)
  }
  return(DT_subset)
}
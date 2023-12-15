desc = function(x) {
  
  A = c(mean(x, na.rm = T), median(x, na.rm = T), sd(x, na.rm = T),
        quantile(x, c(0.01, 0.05, 0.25, 0.75, 0.95, 0.99), na.rm = T))
  names(A) = c('Mean', 'Median', 'S.D.', '1st', '5th', '25th', '75th', '95th', '99th')
  round(A, 2)
  
}

descB = function(x) {
  
  A = c(mean(x, na.rm = T), median(x, na.rm = T), min(x, na.rm = T), max(x, na.rm = T), length(x))
  names(A) = c('Mean', 'Med.', 'Min.', 'Max.', 'N.Obs.')
  c(formatC(A[1:4], digits = 2, format = 'f'), formatC(A[5], digits = 0, format = 'f'))
  
}

descfull = function(x, dig = 2) {
  
  A = c(mean(x, na.rm = T), median(x, na.rm = T), sd(x, na.rm = T), min(x, na.rm = T), 
        quantile(x, c(0.05, 0.25, 0.75, 0.95), na.rm = T), max(x, na.rm = T),
        length(x[!is.na(x)]))
  names(A) = c('Mean', 'Median', 'S.D.', 'Min', '5th', '25th', '75th', '95th', 'Max', 'N.Obs.')
  c(formatC(A[1:9], digits = dig, format = 'f'), formatC(A[10], digits = 0, format = 'f'))
  
}

# define winsorization function
winsorize = function(x, lower, upper) {
  
  q_low = quantile(x, lower, na.rm = T)
  q_upp = quantile(x, upper, na.rm = T)
  x[x<q_low & !is.na(x)] = q_low
  x[x>q_upp & !is.na(x)] = q_upp
  x
  
}

# define truncation function
truncate = function(x, lower, upper) {
  
  q_low = quantile(x, lower, na.rm = T)
  q_upp = quantile(x, upper, na.rm = T)
  x[x<q_low & !is.na(x)] = NA
  x[x>q_upp & !is.na(x)] = NA
  x
  
}
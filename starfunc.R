starfunc = function(x, p, th) {
  
  ast = rep('~~~', length(x))
  ast[p <= 0.10] = '*~~'
  ast[p <= 0.05] = '**~'
  ast[p <= 0.01] = '***'
  
  A = paste0('~', formatC(x, digits = th, format = 'f'), ast)
  A[x >= 0] = paste0('~', A[x >= 0])
  A
  
}

sefunc = function(x, th) {
  
  A = paste0('(', formatC(x, digits = th, format = 'f'), ')~~')
  A[x >= 0] = paste0('~', A[x >= 0])
  A
  
}
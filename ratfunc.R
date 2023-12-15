# define unified rating fucntion
rat = function(x) {
  if (x$RATING %in% c('AAA', 'Aaa')) {'AAA'} 
  else if (x$RATING %in% c('AA+', 'AA', 'AA-', 'Aa3', 'Aa1' , 'Aa2', 'Aa')) {'AA'}
  else if (x$RATING %in% c('A', 'A-', 'A+', 'A1', 'A2', 'A3')) {'A'}
  else if (x$RATING %in% c('BBB-', 'BBB', 'BBB+', 'Baa', 'Baa1', 'Baa2', 'Baa3')) {'BBB'}
  else if (x$RATING %in% c('BB', 'BB-', 'BB+', 'Ba1', 'Ba2', 'Ba3')) {'BB'}
  else if (x$RATING %in% c('B', 'B-', 'B+', 'B1', 'B2', 'B3')) {'B'}
  else if (x$RATING %in% c('CCC', 'CCC-', 'CCC+', 'Caa1', 'Caa2', 'Caa3')) {'CCC'}
  else if (x$RATING %in% c('CC', 'Ca')) {'CC'}
  else if (x$RATING %in% c('C')) {'C'}
  else if (x$RATING %in% c('DD', 'D')) {'D'}
  else {'NR'}
}

# define unified numeric rating fucntion
ratnum = function(x) {
  if (x$RATING %in% c('AAA', 'Aaa')) {1}
  else if (x$RATING %in% c('AA+','Aa1')) {2}
  else if (x$RATING %in% c('AA' ,'Aa2')) {3}
  else if (x$RATING %in% c('AA-','Aa3')) {4}
  else if (x$RATING %in% c('A+','A1')) {5}
  else if (x$RATING %in% c('A' ,'A2')) {6}
  else if (x$RATING %in% c('A-','A3')) {7}
  else if (x$RATING %in% c('BBB+','Baa1')) {8}
  else if (x$RATING %in% c('BBB' ,'Baa', 'Baa2')) {9}
  else if (x$RATING %in% c('BBB-','Baa3')) {10}
  else if (x$RATING %in% c('BB+','Ba1')) {11}
  else if (x$RATING %in% c('BB' ,'Ba2')) {12}
  else if (x$RATING %in% c('BB-','Ba3')) {13}
  else if (x$RATING %in% c('B+','B1')) {14}
  else if (x$RATING %in% c('B' ,'B2')) {15}
  else if (x$RATING %in% c('B-','B3')) {16}
  else if (x$RATING %in% c('CCC+','Caa1')) {17}
  else if (x$RATING %in% c('CCC' ,'Caa2')) {18}
  else if (x$RATING %in% c('CCC-','Caa3')) {19}
  else if (x$RATING %in% c('CC', 'Ca')) {20}
  else if (x$RATING %in% c('C')) {21}  
  else if (x$RATING %in% c('DD', 'D')) {22}
  else {23}
}
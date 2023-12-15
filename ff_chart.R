# ==================================================================================================
# [] PLOT A FF CHART
#  - RUN TIME: 
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven); library(readxl)
library(xtable); library(stargazer); library(lfe); library(plm); library(extrafont)

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

# read data
ff = read_excel('./ff/ff_clean.xlsx', sheet = 2)
ff = as.data.table(ff)

# plot for text
namef = '../Apps/ShareLaTeX/Bond reversals/ff.pdf'

pdf(file = namef, family = 'CM Roman', width = 16, height = 6, pointsize = 26)

par(mar = c(2,2,0.5,0.5), xpd = F)

plot(x = ff[year >= 2007, year], y = ff[year >= 2007, ratiosmart*100], type = 'n',
     las = 1, xlab = '', ylab = '', xaxt = 'n', ylim = c(10, 40))
axis(1, at = 2007:2018, labels = 2007:2018, tick = T)

abline(h = seq(10, 40, 5), lwd = 0.5, lty = 2, col = 'lightgray')
abline(v = seq(2007, 2018, 1), lwd = 0.5, lty = 2, col = 'lightgray')
lines(x = ff[year >= 2007, year], y = ff[year >= 2007, ratiosmart*100], type = 'b', pch = 19,
      col = 'black', lwd = 3)

text(x = 2007, y = 39, labels = '%', adj = 1)

dev.off()
embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))

# plot for slides
namef = '../Apps/ShareLaTeX/Bond reversals/ff_slides.pdf'

pdf(file = namef, family = 'CM Sans', width = 16, height = 6, pointsize = 26)

par(mar = c(2,2,0.5,0.5), xpd = F)

plot(x = ff[year >= 2007, year], y = ff[year >= 2007, ratiosmart*100], type = 'n',
     las = 1, xlab = '', ylab = '', xaxt = 'n', ylim = c(10, 40))
axis(1, at = 2007:2018, labels = 2007:2018, tick = T)

abline(h = seq(10, 40, 5), lwd = 0.5, lty = 2, col = 'lightgray')
abline(v = seq(2007, 2018, 1), lwd = 0.5, lty = 2, col = 'lightgray')
lines(x = ff[year >= 2007, year], y = ff[year >= 2007, ratiosmart*100], type = 'b', pch = 19,
      col = 'black', lwd = 3)

text(x = 2007, y = 39, labels = '%', adj = 1)

dev.off()
embed_fonts(namef, outfile=gsub('.pdf', '_emb.pdf', sub('./fig/', './slides/', namef)))
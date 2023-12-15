# ==================================================================================================
# Update of the stylized chart
# ==================================================================================================

# load data.table package
library(data.table); library(magrittr); library(fastmatch)
library(pbapply); library(xts); library(nleqslv); library(doParallel); library(iterators)
library(timeDate); library(bizdays); library(fasttime); library(haven)
library(xtable); library(stargazer); library(lfe); library(plm); library(extrafont)
library(sandwich); library(lmtest); library(car); library(plotrix)
library(Hmisc); library(estimatr); library(remotes)
library(cowplot); library(grid); library(ggridges); library(scales)
library(colorspace)

# color-blind palette
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
         "#888888")
bdf_pale <- rgb(0.192, 0.549, 0.906, alpha = 0.35)
bdf <- rgb(0.192, 0.549, 0.906, alpha = 1)

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.55.0/bin/gswin64c.exe")

# set locale
Sys.setlocale(locale = 'US')

# load delete function
source('./codes/delete.R')
source('./codes/desc.R')
source('./codes/ratfunc.R')
source('./codes/corstars.R')

# path
namef = '../Apps/ShareLaTeX/Bond reversals/stylized_2023.pdf'

# Data Table for TS
dt1 <- data.table(x = as.numeric(0:60))
dt1[x %in% 0:30, Low := x]
dt1[x %in% 31:60, Low := 40-x/3]
dt1[x %in% 31:60, High := 35-x/6]
dt1[x == 30, High := Low]
dt1 <- melt(dt1, id.vars = 'x')

p1 <- dt1 %>%
  ggplot(mapping = aes(x = x, y = value)) +
  ggtitle(label = 'Individual bond price path conditional on trading volume') +
  geom_line(aes(linetype = variable), lwd = 1.25) + 
  geom_point(aes(x=0,y=0),colour="black", shape = 21, fill = 'white', stroke = 2) + 
  geom_point(aes(x=30,y=30),colour="black", shape = 21, fill = 'white', stroke = 2) + 
  geom_point(aes(x=60,y=20),colour="black", shape = 21, fill = 'white', stroke = 2) + 
  geom_point(aes(x=60,y=25),colour="black", shape = 21, fill = 'white', stroke = 2) + 
  scale_x_continuous(name = element_blank(), breaks = c(0, 30, 60),
                     labels = c('Day 0', 'Day 1', 'Day 2'), limits = c(0,60)) +
  scale_y_continuous(name = element_blank(), breaks = c(0, 20, 25, 30),
                     labels = c(expression(P[0]), expression(P[2]^{Low}), expression(P[2]^{High}),
                                expression(P[1])),
                     limits = c(0,30)) +
  labs(linetype ='Day 1 volume:') +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "lightgray"), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.18), legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  annotate("segment", x = 60, xend = 60, y = 20.5, yend = 24.5, colour = "darkgray",
           arrow = arrow(ends = "both", angle = 30, length = unit(.2,"cm")), size=2) +
  annotate(geom = "curve", x = 50, y = 18, xend = 60, yend = 22.5,
    curvature = .2, arrow = arrow(length = unit(0.2, "cm")), size=2) +
  annotate(geom = "text", x = 40, y = 18, label = "Volume offset", hjust = "left", size = 4)

p1

# data table for CS
dt2 <- data.table(x = as.numeric(c(0,30,60)))
dt2[, Dealer := 15 - x/12]
dt2[, Client := 14 + x/6]
dt2 <- melt(dt2, id.vars = 'x')

p2 <- dt2 %>%
  ggplot(mapping = aes(x = x, y = value, group = variable)) +
  ggtitle(label = 'Volume offset in the cross-section of bonds with differnt info asymmetry') +
  geom_line(color = "black", lwd = 1.25) +
  geom_point(aes(shape = variable), size = 3, color = 'black', stroke = 2, fill = 'white') + 
  scale_shape_manual(values = c(22, 24)) +
  #geom_point(aes(x=0,y=0),colour="black", shape = 21, fill = 'white', stroke = 2) + 
  scale_x_continuous(name = element_blank(), breaks = c(0, 30, 60),
                     labels = c('Low', 'Median', 'High'), limits = c(0,60)) +
  scale_y_continuous(name = element_blank(), breaks = c(14, 19, 24),
                     labels = c(expression(0.75~"M"), expression("M"), expression(1.25~"M")),
                     limits = c(0,30)) +
  labs(shape ='Who provides liquidity:') +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "lightgray"), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.18), legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
p2


p <- plot_grid(p1, p2, nrow = 1)
p

ggsave(plot = p, filename = namef, width = 16, height = 6)
library(tidyverse)
library(gridExtra)
g1<-ggplot(mpg, aes(displ, hwy)) + geom_point() + 
  geom_smooth(span=1.0)+
  theme_bw() + theme(panel.grid=element_blank())

g2<-ggplot(mpg, aes(displ, hwy)) + geom_point() + 
  geom_smooth(span=0.75)+
  theme_bw() + theme(panel.grid=element_blank())

g3<-ggplot(mpg, aes(displ, hwy)) + geom_point() + 
  geom_smooth(span=0.50)+
  theme_bw() + theme(panel.grid=element_blank())

g4<-ggplot(mpg, aes(displ, hwy)) + geom_point() + 
  geom_smooth(span=0.25)+
  theme_bw() + theme(panel.grid=element_blank())

grid.arrange(g1, g2, g3, g4, ncol=2)

## Written by: US EPA, National Center for Environmental Economics; July 2021

##########################
#################  LIBRARY
##########################

## Clear worksace
rm(list = ls())
gc()

## This function will check if a package is installed, and if not, install it
pkgTest <- function(x) {
  if (!require(x, character.only = T)) {
    install.packages(x, dep = T)
    if(!require(x, character.only = T)) stop("Package not found")
  }
}

## These lines load the required packages
packages <- c('magrittr','tidyverse','ggplot2','ggrepel','wesanderson','stringr','readxl') ## you can add more packages here
lapply(packages, pkgTest)

options(scipen=999) # remove scientific notation

##########################
#################### Build 
##########################

page <- read_excel('page_gdp_loss.xlsx') %>% group_by(year) %>% 
  summarise(temp=mean(temp),gdp=sum(gdp),damages=sum(damages),gdp_loss=mean(gdp_loss_perc)) %>% 
  filter(temp<5) %>% 
  mutate(labels=case_when(temp==max(temp)~'Kikstra et al. 2021 \nPAGE2020',TRUE~''))

## TABLE 2 COLUMN 4 LESS THAN 4 DEG
## t2 |       .5953827    .1904833      3.13   0.014     .1561274    1.034638
## mkt_t2 |  -.6217335    .2261748     -2.75   0.025    -1.143293   -.1001736
## cat_t2 |   .2598511    .2672706      0.97   0.359     -.356476    .8761783
## prod_t2 |  .1133249    .125483       0.90   0.393    -.1760394    .4026892
## cross |     1.7005     .3426846      4.96   0.001     .9102677    2.490732

## TABLE 2 COLUMN 8 GREATER THAN 4 DEG
## t2 |         .3181497   .1021783     3.11   0.011     .0904823    .5458172
## mkt_t2 |    -.3445005   .1558482    -2.21   0.052    -.6917519    .0027508
## cat_t2 |     .3622743   .1031227     3.51   0.006     .1325026     .592046
## prod_t2 |    .3982305   .2373944     1.68   0.124    -.1307171    .9271781
## cross |     1.7005      .3306892     5.14   0.000     .9636784    2.437321

## temp vector
temps = tibble(temp=seq(0,4.7,0.1))

# ## damage functions

# ## TEMPS FOR HOWARD AND STERNER LESS THAN 4 DEG
# dams = temps %>%
#   mutate(`Tol (2009)`        = 0.00267*temp^2,
#           `Nordhaus and Moffat (2017) \nDICE2016`        = 0.00236*temp^2,
#           `Nordhaus and Boyer (2009) \nDICE2010`        = 0.00205*temp^2,
#           # `Nordhaus (2019)`        = (0.00236*3.5)*temp^2,
#           # `Weitzman (2012)`                              = (0.00236*temp^2 + 5.07e-6*temp^6.754)/(1+0.00236*temp^2 + 5.07e-6*temp^6.754),
#           # `H&S (2017) \nMarket only`      = (0.5953827/100)*temp^2,
#           `H&S (2017) \nNon-catastrophic` = ((0.5953827*1.25)/100)*temp^2,
#           `H&S  (2017) \n+ Catastrophic` = ((0.5953827*1.25+0.2598511)/100)*temp^2,
#           `H&S  (2017) (2017) \n+ Productivity`     = (((0.5953827+0.1133249)*1.25)/100)*temp^2)

# # # # ## ALL TEMPS FOR HOWARD AND STERNER, INCLUDING GREATER THAN 4 DEG
dams = temps %>%
  mutate(`Tol (2009)`        = 0.00267*temp^2,
         `Nordhaus and Moffat (2017) \nDICE2016`        = 0.00236*temp^2,
         `Nordhaus and Boyer (2009) \nDICE2010`        = 0.00205*temp^2,
         # `Nordhaus (2019)`        = (0.00236*3.5)*temp^2,
         # `Weitzman (2012)`                              = (0.00236*temp^2 + 5.07e-6*temp^6.754)/(1+0.00236*temp^2 + 5.07e-6*temp^6.754),
         # `Howard and Sterner (2017) \nMarket only`      = (0.3181497/100)*temp^2,
         `H&S (2017) \nNon-catastrophic` = ((0.3181497*1.25)/100)*temp^2,
         `H&S  (2017) \n+ Catastrophic`     = (((0.3181497)*1.25+0.3622743)/100)*temp^2,
         `H&S  (2017) \n+ Productivity`     = (((0.3181497+0.3982305)*1.25)/100)*temp^2)


dams %<>% pivot_longer(., cols=!temp, names_to='damage_fun', values_to='gdp_loss')

## add uncertainty
## Howard and Sterner, Table 2, Column 4, preferred less than 4 deg
#     e(V) |         t2      mkt_t2      cat_t2     prod_t2       cross 
# -------------+------------------------------------------------------------
#   t2 |       0.03628389                                                 
#   mkt_t2 |  -0.03628389   0.05115502                                     
#   cat_t2 |  -0.04206286   0.04206286   0.07143358                         
#   prod_t2 | -1.274e-17   -0.01407326   2.389e-17   0.01574598             
#   cross |   -2.127e-17   -0.01413532   4.287e-17   0.00436625   0.11743273 

## Howard and Sterner, Table 2, Column 4, preferred less than 4 deg
# e(V) |         t2             mkt_t2      cat_t2     prod_t2       cross 
# -------------+------------------------------------------------------------
#   t2 |         0.01044041                                                 
#   mkt_t2 |    -0.01044041   0.02428865                                     
#   cat_t2 |    -0.01047679   0.01047679   0.01063429                         
#   prod_t2 |    7.578e-18   -0.01384824  -1.015e-17   0.05635608             
#   cross |      1.530e-17   -0.01316305  -1.831e-17   0.01316305   0.10935535 

dams_unc = temps %>% 
  mutate(`DICE2016`                                     = 0.00118*temp^2)
dams_unc %<>% pivot_longer(., cols=!temp, names_to='damage_fun', values_to='gdp_loss_sd')


data = left_join(dams,dams_unc) %>% mutate(gdp_loss_sd=case_when(is.na(gdp_loss_sd)~0,TRUE~gdp_loss_sd))

# a1 = 0.0127
# a2 = -0.0005
# temp = get_scc_path() %>%
#   mutate(Damages="Burke et al (2015) - PAGE-ICE/2020")
# results = rbind(results,temp)
# 
# a1 = -0.001126
# a2 = 0.000818
# temp = get_scc_path() %>%
#   mutate(Damages="Waldhoff et al. (2014) - FUND3.9")
# results = rbind(results,temp)

##########################
##################### Plot 
##########################

data %<>% mutate(labels=case_when(temp==4.7~damage_fun,TRUE~''),
                 colors=str_extract(damage_fun, '[A-Za-z]+'))

data %>% 
  ggplot() +
  geom_line(aes(x=temp,y=gdp_loss,group=damage_fun,color=damage_fun)) +
  geom_line(data=page,aes(x=temp,y=gdp_loss)) +
  # geom_ribbon(aes(x=temp, ymin=gdp_loss-(1.64485*gdp_loss_sd), ymax=gdp_loss+(1.64485*gdp_loss_sd), fill=colors, group=damage_fun), color=NA, linetype='dotted',alpha = 0.05) +
  geom_label_repel(aes(x=temp,y=gdp_loss,color=damage_fun,label=labels), size=4,max.overlaps=100,nudge_x =0.5) +
  geom_label_repel(data=page,aes(x=temp,y=gdp_loss,label=labels), size=4,max.overlaps=100,nudge_x =0.5) +
  scale_color_manual(values=wes_palette(name="BottleRocket1")) +
  scale_fill_manual(values=wes_palette(name="BottleRocket1")) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5), limits=c(0,6)) +
  scale_y_continuous(breaks=c(0,0.05,0.10,0.15,0.20,0.25,0.30,0.35),labels=scales::percent_format()) +
  labs(x="Change in Temperature",y="GDP Loss") +
  theme_minimal() + 
  theme(legend.position="none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave("figures/gdp_loss_under_all_temps.png",width=11,height=8)


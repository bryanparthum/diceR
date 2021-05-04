# time periods (5 years per period)
time_horizon = 60
  
# availability of fossil fuels
fosslim  = 6000      # maximum cumulative extraction fossil fuels (gtc)

# time step
tstep    = 5         # years per period
  
# preferences
elasmu   = 1.45      # elasticity of marginal utility of consumption
prstp    = .015      # initial rate of social time preference per year
  
# population and technology
gama     = .300      # capital elasticity in production function
pop0     = 7403      # initial world population 2015 (millions)
popadj   = 0.134     # growth rate to calibrate to 2050 pop projection
popasym  = 11500     # asymptotic population (millions)
dk       = .100      # depreciation rate on capital (per year)
q0       = 105.5     # initial world gross output 2015 (trill 2010 usd)
k0       = 223       # initial capital value 2015 (trill 2010 usd)
a0       = 5.115     # initial level of total factor productivity
ga0      = 0.076     # initial growth rate for tfp per 5 years 
dela     = 0.005     # decline rate of tfp per 5 years
s        = 0.2       # savings rate

# emissions parameters
gsigma1  = -0.0152   # initial growth of sigma (per year)
dsig     = -0.001    # decline rate of decarbonization (per period)
eland0   = 2.6       # carbon emissions from land 2015 (gtco2 per year)
deland   = .115      # decline rate of land emissions (per period)
e0       = 35.85     # industrial emissions 2015 (gtco2 per year)
miu0     = .03       # initial emissions control rate for base case 2015

# carbon cycle
# initial conditions
mat0     = 851       # initial concentration in atmosphere 2015 (gtc)
mu0      = 460       # initial concentration in upper strata 2015 (gtc)
ml0      = 1740      # initial concentration in lower strata 2015 (gtc)
mateq    = 588       # equilibrium concentration atmosphere  (gtc)
mueq     = 360       # equilibrium concentration in upper strata (gtc)
mleq     = 1720      # equilibrium concentration in lower strata (gtc)

# flow paramaters
b12      = .12       # carbon cycle transition matrix
b23      = 0.007     # carbon cycle transition matrix
  
# these are for declaration and are defined later
# b11      carbon cycle transition matrix
# b21      carbon cycle transition matrix
# b22      carbon cycle transition matrix
# b32      carbon cycle transition matrix
# b33      carbon cycle transition matrix
# sig0     carbon intensity 2010 (kgco2 per output 2005 usd 2010)

# climate model parameters
t2xco2   = 3.1       # equilibrium temp impact (oc per doubling co2)
fex0     = 0.5       # 2015 forcings of non-co2 ghg (wm-2)
fex1     = 1.0       # 2100 forcings of non-co2 ghg (wm-2)
tocean0  = .0068     # initial lower stratum temp change (c from 1900)
tatm0    = 0.85      # initial atmospheric temp change (c from 1900)
c1       = 0.1005    # climate equation coefficient for upper level
c3       = 0.088     # transfer coefficient upper to lower stratum
c4       = 0.025     # transfer coefficient for lower level
fco22x   = 3.6813    # forcings of equilibrium co2 doubling (wm-2)
  
# climate damage parameters
a10      = 0         # initial damage intercept
# a20                # initial damage quadratic term
a1       = 0         # damage intercept
a2       = 0.00236   # damage quadratic term
a3       = 2.00      # damage exponent
a4       = 5.07e-6   # Weitzman damage term
a5       = 6.754     # damage exponent Weitzman 2012

# abatement cost
expcost2 = 2.6       # exponent of control cost function
pback    = 550       # cost of backstop 2010$ per tco2 2015
gback    = .025      # initial cost decline backstop cost per period
limmiu   = 1.2       # upper limit on control rate after 2150
tnopol   = 45        # period before which no emissions controls base
cprice0  = 2         # initial base carbon price (2010$ per tco2)
gcprice  = .02       # growth rate of base carbon price per year
  
# scaling and inessential parameters
# note that these are unnecessary for the calculations
# they ensure that mu of first period's consumption =1 and pv cons = pv utilty
scale1   = 0.0302455265681763   # multiplicative scaling coefficient
scale2   = -10993.704     # additive scaling coefficient

# parameters
# l[t]          level of population and labor
# al[t]         level of total factor productivity
# sigma[t]      co2-equivalent-emissions output ratio
# rr[t]         average utility social discount rate
# ga[t]         growth rate of productivity from
# forcoth[t]    exogenous forcing for other greenhouse gases
# gl[t]         growth rate of labor
# gcost1        growth of cost factor
# gsig[t]       change in sigma (cumulative improvement of energy efficiency)
# etree[t]      emissions from deforestation
# cumetree[t]   cumulative from land
# cost1[t]      adjusted cost for backstop
# lam           climate model parameter
# gfacpop[t]    growth factor population
# pbacktime[t]  backstop price
# optlrsav      optimal long-run savings rate used for transversality
# scc[t]        social cost of carbon
# cpricebase[t] carbon price in base case
# photel[t]     carbon price under no damages (hotelling rent condition)
# ppm[t]        atmospheric concentrations parts per million
# atfrac[t]     atmospheric share since 1850
# atfrac2010[t]     atmospheric share since 2010 ;

# parameters for long-run consistency of carbon cycle
b11 = 1 - b12
b21 = b12*mateq/mueq
b22 = 1 - b21 - b23
b32 = b23*mueq/mleq
b33 = 1 - b32 

# further definitions of parameters
a20 = a2
sig0 = e0/(q0*(1-miu0))
lam = fco22x/t2xco2

#
ga = ga0*exp(-dela*5*(0:(time_horizon-1)))

#
etree = eland0*(1-deland)^(0:(time_horizon-1))

#
forcoth = array(fex0+(1/17)*(fex1-fex0),time_horizon)
forcoth[-(1:17)] = (fex1-fex0)

pbacktime = pback*(1-gback)^(0:(time_horizon-1))

l = array(pop0,time_horizon)
al = array(a0,time_horizon) 
gsig = array(gsigma1,time_horizon)
sigma = array(sig0,time_horizon)
cumetree = array(100,time_horizon)
cost1 = array(pback*sig0/expcost2/1000,time_horizon)

for (t in 2:time_horizon) {
  
  l[t] = l[t-1]*(popasym/l[t-1])^popadj
  
  al[t] = al[t-1]/(1-ga[t-1])
  
  gsig[t] = gsig[t-1]*(1+dsig)^tstep
  
  sigma[t] = sigma[t-1]*exp(gsig[t-1]*tstep)
  
  cost1[t] = pbacktime[t]*sigma[t]/expcost2/1000

}

run_dice = function(perturbation_year=-1,damfun) {

  # variables
  # damfrac[t]      damages as fraction of gross output
  # ygross[t]       gross world product gross of abatement and damages (trillions 2005 usd per year)
  # ynet[t]         output net of damages equation (trillions 2005 usd per year)
  # abatecost[t]    cost of emissions reductions  (trillions 2005 usd per year)
  # y[t]            gross world product net of abatement and damages (trillions 2005 usd per year)
  # i[t]            investment (trillions 2005 usd per year)
  # c[t]            consumption (trillions 2005 us dollars per year)
  # k[t]            capital stock (trillions 2005 us dollars)
  # eind[t]         industrial emissions (gtco2 per year)
  # e[t]            total co2 emissions (gtco2 per year)
  # cca[t]          cumulative industrial carbon emissions (gtc)
  # forc[t]         increase in radiative forcing (watts per m2 from 1900)
  # mat[t]          carbon concentration increase in atmosphere (gtc from 1750)
  # ml[t]           carbon concentration increase in lower oceans (gtc from 1750)
  # mu[t]           carbon concentration increase in shallow oceans (gtc from 1750)
  # tatm[t]         increase temperature of atmosphere (degrees c from 1900)
  # tocean[t]       increase temperatureof lower oceans (degrees c from 1900)
  # miu[t]          emission control rate ghgs
  
  # intialize space for state variables and set starting values
  damfrac   = array(NA,      time_horizon)
  ygross    = array(NA,      time_horizon)
  ynet      = array(NA,      time_horizon)
  abatecost = array(NA,      time_horizon)
  y         = array(NA,      time_horizon)
  i         = array(NA,      time_horizon)
  c         = array(NA,      time_horizon)
  k         = array(k0,      time_horizon)
  eind      = array(NA,      time_horizon)
  e         = array(NA,      time_horizon)
  cca       = array(400,     time_horizon)
  forc      = array(NA,      time_horizon)
  mat       = array(mat0,    time_horizon)
  ml        = array(ml0,     time_horizon)
  mu        = array(mu0,     time_horizon)
  tatm      = array(tatm0,   time_horizon)
  tocean    = array(tocean0, time_horizon)
  miu       = array(miu0,    time_horizon)
  
  #
  for (t in 1:time_horizon) {
    
    # equation for damage fraction
    if (missing(damfun)) {
      damfrac[t] = a1*tatm[t]+a2*tatm[t]^a3
      dam_fun = "DICE2016"
    }
    else if (damfun=="Weitzman") {
      damfrac[t] = a1*tatm[t]+a2*tatm[t]^a3+a4*tatm[t]^a5
      dam_fun = "Weitzman (2012)"
    }
    
    # output gross equation
    ygross[t] = (al[t]*(l[t]/1000)^(1-gama))*(k[t]^gama)
    
    # output net of damages equation
    ynet[t] = ygross[t]*(1-damfrac[t])
    
    # cost of emissions reductions equation
    abatecost[t] = ygross[t]*cost1[t]*miu[t]^expcost2
    
    # output net equation
    y[t] = ynet[t]-abatecost[t]
    
    # savings rate equation
    i[t] = s*y[t]
    
    # consumption equation
    c[t] = y[t]-i[t]
    
    # capital balance equation
    k[t+1] = (1-dk)^tstep*k[t]+tstep*i[t]
    
    # industrial emissions
    eind[t] = sigma[t]*ygross[t]*(1-(miu[t]))
    
    if (perturbation_year==(2015+t*tstep))
      eind[t] = eind[t]+1
    
    # emissions equation
    e[t] = eind[t]+etree[t]
    
    # cumulative industrial carbon emissions
    cca[t+1] = cca[t]+eind[t]*5/3.666
    
    # atmospheric concentration equation
    mat[t+1] = mat[t]*b11+mu[t]*b21+e[t]*(5/3.666)
    
    # lower ocean concentration
    ml[t+1] = ml[t]*b33+mu[t]*b23
    
    # shallow ocean concentration
    mu[t+1] = mat[t]*b12+mu[t]*b22+ml[t]*b32
    
    # radiative forcing equation
    forc[t+1] = fco22x*((log((mat[t+1]/588.000))/log(2)))+forcoth[t+1]
    
    # temperature-climate equation for atmosphere
    tatm[t+1] = tatm[t]+c1*(forc[t+1]-(fco22x/t2xco2)*tatm[t]-c3*(tatm[t]-tocean[t]))

    # temperature-climate equation for lower oceans
    tocean[t+1] = tocean[t]+c4*(tatm[t]-tocean[t])
    
  }

  # return(c)
  return(data.frame("consumption"= c,'abatecost'=abatecost,'tatm'=tatm[1:60],'emissions'=e,'ygross'=ygross,'ynet'=ynet, 'capital'=k[1:60], 'savings'=i[1:60], 'dam_frac'=damfrac))
}

library(tidyverse)
library(magrittr)

base = run_dice()
pert = run_dice(perturbation_year=2030)
diff = pert-base 
diff %<>% mutate(dam_fun = 'DICE2016',year=seq(2015,2310,5))

base_w = run_dice(damfun="Weitzman")
pert_w = run_dice(perturbation_year=2030,damfun="Weitzman")
diff_w = pert_w-base_w 
diff_w %<>% mutate(dam_fun = 'Weitzman (2012)',year=seq(2015,2310,5))

bases <- rbind(base %>% mutate(dam_fun = 'DICE2016',year=seq(2015,2310,5)),base_w %>% mutate(dam_fun = 'Weitzman (2012)',year=seq(2015,2310,5)))
diffs <- rbind(diff,diff_w)

library(ggplot2)
library(gridExtra)
library(cowplot)
options(scipen=999)

d <- ggplot(diffs) + 
  geom_line(aes(y=dam_frac,x=year,color=dam_fun)) +
  labs(title="Damage Fraction",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

s <- ggplot(diffs) + 
  geom_line(aes(y=savings,x=year,color=dam_fun)) +
  labs(title="Savings",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

k <- ggplot(diffs) + 
  geom_line(aes(y=capital,x=year,color=dam_fun)) +
  labs(title="Capital",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

c <- ggplot(diffs) + 
  geom_line(aes(y=consumption,x=year,color=dam_fun)) +
  labs(title="Consumption",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

a <- ggplot(diffs) + 
  geom_line(aes(y=abatecost,x=year,color=dam_fun)) +
  labs(title="Abatement Cost",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")


t <- ggplot(diffs) + 
  geom_line(aes(y=tatm,x=year,color=dam_fun)) +
  labs(title="Temperature",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

e <- ggplot(diffs) + 
  geom_line(aes(y=emissions,x=year,color=dam_fun)) +
  labs(title="Emissions",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")


yg <- ggplot(diffs) + 
  geom_line(aes(y=ygross,x=year,color=dam_fun)) +
  labs(title="Y-gross",
       y='',
       x='Year') +
  theme_minimal() + theme(legend.position = "none")


yn <- ggplot(diffs) + 
  geom_line(aes(y=ynet,x=year,color=dam_fun)) +
  labs(title="Y-net",
       y='',
       x='Year') +
  theme_minimal() + theme(legend.position = "none")

legend <- get_legend(
  ggplot(diffs) + 
  geom_line(aes(y=ynet,x=year,color=dam_fun)) +
  labs(color="Difference From CO[2] Pulse in 2030 (1MT):") + theme(legend.position='bottom',legend.key = element_rect(color = NA, fill = NA)))

b <- ggplot()+theme_void()

plot_grid(d,s,k,c,a,t,e,yg,yn,b,legend,b, nrow=4,rel_heights = c(3/10, 3/10, 3/10, 1/10))
ggsave("differences.png",width=8,height=6)


## BASELINES
d <- ggplot(bases) + 
  geom_line(aes(y=dam_frac,x=year,color=dam_fun)) +
  labs(title="Damage Fraction",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

s <- ggplot(bases) + 
  geom_line(aes(y=savings,x=year,color=dam_fun)) +
  labs(title="Savings",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

k <- ggplot(bases) + 
  geom_line(aes(y=capital,x=year,color=dam_fun)) +
  labs(title="Capital",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

c <- ggplot(bases) + 
  geom_line(aes(y=consumption,x=year,color=dam_fun)) +
  labs(title="Consumption",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

a <- ggplot(bases) + 
  geom_line(aes(y=abatecost,x=year,color=dam_fun)) +
  labs(title="Abatement Cost",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")


t <- ggplot(bases) + 
  geom_line(aes(y=tatm,x=year,color=dam_fun)) +
  labs(title="Temperature",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")

e <- ggplot(bases) + 
  geom_line(aes(y=emissions,x=year,color=dam_fun)) +
  labs(title="Emissions",
       y='',
       x='') +
  theme_minimal() + theme(legend.position = "none")


yg <- ggplot(bases) + 
  geom_line(aes(y=ygross,x=year,color=dam_fun)) +
  labs(title="Y-gross",
       y='',
       x='Year') +
  theme_minimal() + theme(legend.position = "none")


yn <- ggplot(bases) + 
  geom_line(aes(y=ynet,x=year,color=dam_fun)) +
  labs(title="Y-net",
       y='',
       x='Year') +
  theme_minimal() + theme(legend.position = "none")

legend <- get_legend(
  ggplot(diffs) + 
    geom_line(aes(y=ynet,x=year,color=dam_fun)) +
    labs(color="Baseline:") + theme(legend.position='bottom',legend.key = element_rect(color = NA, fill = NA)))

b <- ggplot()+theme_void()

plot_grid(d,s,k,c,a,t,e,yg,yn,b,legend,b, nrow=4,rel_heights = c(3/10, 3/10, 3/10, 1/10))
ggsave("baseline.png",width=8,height=6)


# perturbation_year = 2030
# c_perturb = run_dice(perturbation_year)
# 
# c_diff = (c_base-c_perturb)[-(1:((perturbation_year-2015)/tstep+1))]
# 
# print(sum(c_diff*tstep/(1+.03)^(0:(length(c_diff)-1)*tstep))*1e12/1e9/5)


# ## testing damage function
# 
# consumption_weitz <- data.frame(weitz=run_dice(damfun="Weitzman"),
#                                 damage_function='Weitzman (2012)',
#                                 year = seq(2015,2310,5))
# consumption_dice  <- data.frame(dice=run_dice(),
#                                 damage_function='DICE2016',
#                                 year = seq(2015,2310,5))
# dice <- rbind(consumption_weitz,consumption_dice)
# 
# library(ggplot2)
# ggplot(c) + 
#   geom_line(aes(y=consumption,x=year,color=damage_function)) +
#   theme_bw()

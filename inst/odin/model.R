## Neil's dengue model transcribed from BM to odin
## Vaccination removed
##

  PI_C <- 3.14159265359 ## PI
  
  DT  <-  1.0  ## timestep in days - keep at 1
  TIME <- step*DT
  YL <- 365 ## year length in days

  NUM_YEAR_ACCUM <- 4 ## time to accumulate incidence for a couple of output variables
  
  ## Age categories
  
  age_per <- YL ## how frequently aging occurs - needs to be year length for precise match to demography
  N_age <- 25 ## number of age classes
  N_age_p1 <- N_age + 1
  agec[1:15] <- 1 ## first 15 are 1 year wide
  agec[16:18] <- 5 ## nexr 3 are 5 years
  agec[19:25] <- 10 ### then 10 year age bands to 100
  ## lower boundary of each age class. Extra entry gives upper boundary of last age class 
  ageb[1] <- 0
  ageb[2:(N_age_p1)] <- ageb[i-1]+agec[i-1]
  
  N_sim <- user() # simulate 5 million people (initial pop size)
  
  age_removal_d[,] <- user()
  age_rate_d[,] <- user()
  pop_size_d[] <- user()
  births_d[,] <- user()
  life_expec_d[,] <- user()

  EQUILIB_YEARS <- user() ## years to run model before FIRST_YEAR reached

  CALIB_YEAR <- user() ## calendar year used for R0 calibration against demography (for age-dependent exposure)
  FIRST_YEAR <- user() ## first calendar year in demography data
  LAST_YEAR <- user() ## last calendar year in demography data
  max_rel_year <- LAST_YEAR-FIRST_YEAR ## number of years covered by demography data
  

  ## handle leap years, ugh!
  START_YEAR <- FIRST_YEAR-EQUILIB_YEARS
  ## Day of week of 1 Jan of START_YEAR - from https://en.wikipedia.org/wiki/Determination_of_the_day_of_the_week
  DOW_START0 <- (1+5*((START_YEAR-1)%% 4) + 4*((START_YEAR-1)%% 100) +6*((START_YEAR-1)%% 400))%% 7 
  DOW_START <- if(DOW_START0==0) 7 else DOW_START0
  CUR_DOW <- ((floor(DOW_START+TIME)-1) %% 7) + 1
  initial(out_CUR_DOW) <- DOW_START
  update(out_CUR_DOW) <- CUR_DOW

  
  LEAP_YEAR <- START_YEAR %% 4
  YEAR_OFFSET <- if(LEAP_YEAR==0) 0 else (4.0-LEAP_YEAR)
  TIME_OFFSET <- TIME+YEAR_OFFSET*YL/DT
  DAYS_IN_4Y <- YL*4+1
  FOUR_YEARS <- floor((TIME_OFFSET)*DT / DAYS_IN_4Y)
  DAYS_AFTER_4Y <- TIME_OFFSET*DT-FOUR_YEARS*DAYS_IN_4Y
  YEAR <- FOUR_YEARS*4.0 + (if(DAYS_AFTER_4Y<366.0) DAYS_AFTER_4Y/366.0 else 1.0+(DAYS_AFTER_4Y-366.0)/365.0) - YEAR_OFFSET
  CALENDAR_YEAR <- START_YEAR + YEAR
  ORDINAL_YEAR <- floor(CALENDAR_YEAR+0.0001) # 0.0001 to allow for numerical inaccuracy
  CUR_YEAR_LENGTH <- if(ORDINAL_YEAR %% 4==0) 366 else 365
  DAY_OF_YEAR <- floor((CALENDAR_YEAR-ORDINAL_YEAR)*CUR_YEAR_LENGTH+0.25) # 0.25 to allow for numerical inaccuracy
  
  initial(out_year) <- 0
  update(out_year) <- CALENDAR_YEAR
  
  ### Epi data (not actually read into model - these vars control output)
  ## A lot of work converting to ISO weeks!
  
  DATA_START_YEAR <- user() ## calendar year data starts
  DATA_NUM_WEEKS <- user() ## number of weeks in epi data
  DATA_NUM_YEARS <- floor(DATA_NUM_WEEKS/52) ## accurate for DATA_NUM_YEARS up to at least 100
  DATA_SERO_YEAR <- user() ## year of serological survey
  DATA_REPORTING_DELAY <- user() ## reporting (plus incub per etc) delays in days
  DATA_DOW <- ((CUR_DOW-1+DATA_REPORTING_DELAY) %% 7)+1
  DATA_DAY_OF_YEAR0 <- DAY_OF_YEAR + DATA_REPORTING_DELAY
  DATA_YEAR <- if(DATA_DAY_OF_YEAR0>CUR_YEAR_LENGTH) ORDINAL_YEAR+1 else ORDINAL_YEAR
  DATA_DAY_OF_YEAR <- if(DATA_DAY_OF_YEAR0>CUR_YEAR_LENGTH) DATA_DAY_OF_YEAR0-CUR_YEAR_LENGTH else DATA_DAY_OF_YEAR0
  DATA_YEAR_LENGTH <- if(DATA_YEAR %% 4==0) 366 else 365
  LAST_DATA_YEAR_LENGTH <- if((DATA_YEAR-1) %% 4==0) 366 else 365
  
  THURS_DOY_IN_WEEK0 <- 4-DATA_DOW + DATA_DAY_OF_YEAR
  THURS_DOY_IN_WEEK <- if(THURS_DOY_IN_WEEK0>DATA_YEAR_LENGTH) (THURS_DOY_IN_WEEK0-DATA_YEAR_LENGTH) else if(THURS_DOY_IN_WEEK0<=0) (THURS_DOY_IN_WEEK0+LAST_DATA_YEAR_LENGTH) else THURS_DOY_IN_WEEK0
  CUR_ISO_WEEK <- floor((6+THURS_DOY_IN_WEEK)/7)
  DATA_PERIOD_ACTIVE <- if((DATA_YEAR>=DATA_START_YEAR) && ((DATA_YEAR>DATA_START_YEAR)||(DATA_DAY_OF_YEAR>=4)||(CUR_ISO_WEEK==1))) 1 else 0
  initial(DATA_DAY) <- 0
  update(DATA_DAY) <- if(DATA_PERIOD_ACTIVE==0) DATA_DAY else DATA_DAY+DT
  cur_data_day <- floor(DATA_DAY) ## starts at 0
  cur_data_week <- floor(cur_data_day / 7) + 1 ## starts at 1
  initial(out_curisoweek) <- 0
  update(out_curisoweek) <- CUR_ISO_WEEK
  
  ### code to handle climate data
  
  FIRST_CLIM_YEAR <- user() ## calendar year climate data starts
  LAST_CLIM_YEAR <- user() ## calendar year climate data ends
  REPEAT_CLIM_4YEARS <- user() ## how many climate years (*4) to replicated before/after climate data available

  ## temp data from FIRST_CLIM_YEAR up until (but not including) LAST_CLIM_YEAR
  NTP <- user() ## length of climate timeseries in days
  climate_d[,] <- user() ## cols are temperature, rainfall, in one day steps
  
  rel_first_clim_year <- FIRST_CLIM_YEAR-START_YEAR
  ## days after start when climate data starts, accounting for leap years
  rel_first_clim_day <- (rel_first_clim_year*YL + floor(rel_first_clim_year / 4) + (if(LEAP_YEAR < rel_first_clim_year%%4) 1 else 0))
  rel_last_clim_year <- LAST_CLIM_YEAR-START_YEAR
  ## days after start when climate data ends, accounting for leap years
  rel_last_clim_day <- (rel_last_clim_year*YL + floor(rel_last_clim_year / 4) + (if(LEAP_YEAR < rel_last_clim_year%%4) 1 else 0)) - 1
  max_clim_day <- rel_last_clim_day-rel_first_clim_day
  ## current model day relative to start of climate data
  cur_clim_day <- floor(TIME*DT)-rel_first_clim_day
  initial(out_cur_clim_day) <- 0
  update(out_cur_clim_day) <- cur_clim_day
  ## repeat 4 year chunks of climate data for model times before climate data starts
  norm_rel_clim_day <- (cur_clim_day+REPEAT_CLIM_4YEARS*DAYS_IN_4Y*25) %% (REPEAT_CLIM_4YEARS*DAYS_IN_4Y)
  ## same after (but repeat end of climate timeseries)
  norm_rel_clim_end_day <- max_clim_day-REPEAT_CLIM_4YEARS*DAYS_IN_4Y+(cur_clim_day-max_clim_day) %% (REPEAT_CLIM_4YEARS*DAYS_IN_4Y)
  ## now find row of climate dataset corresponding to current time
  clim_row <- if(cur_clim_day<0) norm_rel_clim_day else if(cur_clim_day>=max_clim_day) norm_rel_clim_end_day else cur_clim_day
  ## find temperature and rainfall for current timestep
  # print("clim_row: {clim_row}")
  temperature <- climate_d[as.integer(clim_row+1),1]
  rainfall  <- climate_d[as.integer(clim_row+1),2] *DT ## if DT<1, spread rain equally across day
  rainfall_mean <- sum(climate_d[,2])/NTP  ## needed to normalise rainfall dependent carrying capacity



  #### dynamic demog
  
  
  YEAR_STEP <- 1 ## time between population size estimates given in demography files
  
  ## current model year relative to FIRST_YEAR
  rel_year <- YEAR - EQUILIB_YEARS  
  
  ## row in age_removal demog file corresponding to current year
  year_row <- 1+(if(rel_year<0) 0 else if(rel_year>max_rel_year) max_rel_year else floor(rel_year))
  
  ## year used to calibrate R0 (for age-dependent exposure)
  year_calib <- 1+CALIB_YEAR-FIRST_YEAR
  
  death[1:N_age] <- age_removal_d[as.integer(year_row),1+i]
  ## correction for first year to allow for continuous births
  death[1] <- 0.4691*death[1]*death[1]+1.9686*death[1]
  ## change time units to days
  deathrt[1:N_age] <- death[i]/YL*DT
  
  ## not really the mean age of people in each age class, but biassed below midpoint to allow for deaths
  mean_age[1:N_age] <- (0.75*ageb[i]+0.25*ageb[i+1])
  
  ## proportion of each age class which age to the next each year
  cur_age_rate[1:N_age] <- age_rate_d[as.integer(year_row),1+i]
  ## age rate is only non-zero once per age_per
  agerts <- if(floor(TIME/age_per) == TIME/age_per) age_per/YL else 0
  ## agerts <- DT/YL
  agert[1:N_age] <- cur_age_rate[i]*agerts  ## /agec[i]
  
  ## initial value used for initialising human population 
  N_init_age0[1:N_age] <- pop_size_d[i+1]
  N_init <- sum(N_init_age0)
  pop_scale <- if(N_sim<=0) 1.0 else N_sim/N_init
  N_init_age[1:N_age] <- pop_size_d[i+1]*pop_scale

  ## currently used for initialising mosquito model
  N_eq <- sum(N_init_age)
  
  ## Life expectancy for each age class - used to calculate YLLs
  life_expec[1:N_age] <- life_expec_d[as.integer(year_row),1+i]
  ## life expectancy used for R0 calibration
  init_lifespan <- life_expec_d[as.integer(year_calib),2]
  init_life_expec[1:N_age] <- life_expec_d[as.integer(year_calib),1+i]
  
  ## births happen continuously
  Nb <- births_d[as.integer(year_row),2]
  births <- pop_scale*Nb*DT/YL
  

  ## proportion of infections causing symptomatic disease in unvaccinated
  
  dis_pri[] <- user() ## prop of primary infectious which are symptomatic
  dis_sec[] <- user()
  dis_tert[] <- user()
  dis_quart[] <- user()
  
  ## proportion of infections causing severe disease
  sdis_pri[] <- user() ## prop of primary infectious which are severe/hosp
  sdis_sec[] <- user()
  sdis_tert[] <- user()
  sdis_quart[] <- user()
  
  
  ## mosquito params
  
  omega <- user()  ## density dep power
  
  
  ## adult mos death rate 
  ## assuming lifespan shows quadratic like temp dependence
  
  delta_p <- user()  ## min death rate (at optimal temperature)
  delta0 <- 0.13 ## approx mean death rate for equilib calculations
  delta_max <- 1.0 ## max death rate (set to 1/day to avoid numerical issues in model)
  ## Temperature dependence in delta
  
  # delta_T0 <- user()  ## min temperature
  # delta_Tw <- user()  ## width of temperature window
  # delta_Tm <- delta_T0 + delta_Tw ## max temperature
  
  delta_Tm <- user()  ## max temp
  delta_fT0 <- user() ## min temp as fraction of max temp
  delta_T0 <- delta_Tm*delta_fT0  ## min temp
  
  delta_pTm <- user() ## power on Tm term
  delta_pT0d <- user() ## difference between power on Tm and power on T0
  delta_pT0 <- delta_pTm+delta_pT0d ## power on T0 term
  delta_Tp <- (delta_T0*delta_pTm+delta_Tm*delta_pT0)/(delta_pT0+delta_pTm) ## temp which gives peak delta_hm
  delta_norm <- 1.0/(((delta_Tp-delta_T0)^delta_pT0)*((delta_Tm-delta_Tp)^delta_pTm)) ## multiplier which give temp dep function max of 1
  ## temperature driven delta
  delta_temperature <- if(temperature <= delta_T0 || temperature >= delta_Tm) 0 else delta_norm * ((temperature-delta_T0)^delta_pT0)*((delta_Tm-temperature)^delta_pTm)
  delta <- 1.0/(1.0/delta_max+delta_temperature*(1.0/delta_p-1.0/delta_max))

  
  epsilon <- user()  ## larval maturation
  Rm <- user() ##  mosq rep no.
  
  ##  Mwt = adult wild type female density per person
  ##  this param can also be changed to give required R0 rather than beta_hm
  # Mwt <- user()  # now used to determined R0
  
  sigma <- user() ##  low density limit of larval death rate
  ## fecundity - computed to give required Rm
  gamma <- Rm*delta0*(epsilon+sigma)/epsilon
  
  ## carrying capacity - computed to give required Mwt
  Kc_mean <- Mwt*delta0*( (epsilon*(gamma-delta0)/(delta0*sigma)-1) ^(-1/omega) )/epsilon
  
  ## rainfall driven carrying capacity
  ##
  ## only dependent on rainfall, not humidity
  ## tau_rain is mean period of accummulation
  ## sat_rain=rainfall (multiplier of mean) accum_rain saturates at (if max_rain -> infty)
  ## max_rain=rainfall (multiplier of mean) level for which 
  ## accum_rain reaches a maximum (declines for higher rainfall due to washout)
  ##
  ## For more complex model, see
  ## https://www.sciencedirect.com/science/article/pii/S0304380018302382
  ## TBC
  ##

  # Kc <- Kc_mean
  
  tau_rain <- user()
  sat_rain <- user()
  max_rain <- user()
  
  ## Model currently normalised so sat_rain and max_rain are relative to mean rain (over timeseries)
  ## This may not be the best paramtererisation over multiple locations
  ## Easy to switch to absolute - just remove /rainfall_mean in next line
  ##
  
  rain_norm <- rainfall/rainfall_mean
  
  ## what we would expect accum_rain to equilibriate to if rainfall was constant
  accum_rain_eq <- 1.0/(1.0+1.0/sat_rain+1.0/(max_rain*max_rain))
  
  ## exp to prevent accum_rain going negative for small max_rain and high rainfall
  change_accum_rain <- exp(-(1+rain_norm/sat_rain+(rain_norm*rain_norm)/(max_rain*max_rain))*DT/tau_rain)
  
  initial(accum_rain) <- accum_rain_eq
  ## implicit DT^2 below intended to allow high rainfall, low max_rain scenario to pull down accum_rain fast
  ## division by rainfall_mean makes accum_rain, sat_rain, max_rain relative to mean rainfall
  update(accum_rain) <- (accum_rain+ DT*rain_norm/tau_rain)*change_accum_rain


  Kc <- Kc_mean*accum_rain/accum_rain_eq
  
  
  kappa <- user()  ## biting rate/day
  
  ## mosquito to human transmission prob (multiplier of temp-dep rather than mean)
  Beta_mh_max <- user() 
  ## human to mosquito transmission prob (multiplier of temp-dep rather than mean)
  Beta_hm_max <- user() 
  
  Beta_mh_mean <- 0.86*Beta_mh_max
  Beta_hm_mean <- 0.86*Beta_hm_max
  
  ## Temperature dependence in beta
  # Beta_T0 <- user()  ## min temperature
  # Beta_Tw <- user()  ## width of temperature window
  # Beta_Tm <- Beta_T0 + Beta_Tw ## max temperature
  # 
  Beta_Tm <- user()  ## max temperature
  Beta_fT0 <- user()  ## min temp as fraction of max temp
  Beta_T0 <- Beta_Tm * Beta_fT0 ## min temperature
  
  Beta_pTm <- user() ## power on Tm term
  Beta_pT0d <- user() ## difference between power on Tm and power on T0
  Beta_pT0 <- Beta_pTm+Beta_pT0d ## power on T0 term
  Beta_Tp <- (Beta_T0*Beta_pTm+Beta_Tm*Beta_pT0)/(Beta_pT0+Beta_pTm) ## temp which gives peak Beta_hm
  Beta_norm <- 1.0/(((Beta_Tp-Beta_T0)^Beta_pT0)*((Beta_Tm-Beta_Tp)^Beta_pTm)) ## multiplier which give temp dep function max of 1
  ## temperature driven beta
  Beta_temperature <- if(temperature <= Beta_T0 || temperature >= Beta_Tm) 0 else (Beta_norm * ((temperature-Beta_T0)^Beta_pT0)*((Beta_Tm-temperature)^Beta_pTm))
  Beta_mh <- Beta_mh_max*Beta_temperature ## assume same form for Beta_hm and Beta_mh
  Beta_hm <- Beta_hm_max*Beta_temperature
  
  ## external force of infection on people
  extInf <- user()
  
  ## Facility to vary human exposure by age
  Acrit <- 3
  FOIagescale <- 1.0
  FOIas[1:N_age] <- if(i>Acrit) FOIagescale else 1
  suscinitpop[1:N_age] <- FOIas[i]*init_life_expec[i]
  R0agescale <- sum(init_life_expec)/sum(suscinitpop)
  initial(out_R0agescale) <- 0
  update(out_R0agescale) <- R0agescale
  
  Mwt <- user() ## required Mwt
  
  Rel_R01 <- user() ## R0 scaling for DENV1
  Rel_R02 <- user() ## R0 scaling for DENV2
  Rel_R03 <- user() ## R0 scaling for DENV3
  Rel_R04 <- user() ## R0 scaling for DENV4
  
  Beta_hm_1 <- Beta_hm*Rel_R01
  Beta_hm_2 <- Beta_hm*Rel_R02
  Beta_hm_3 <- Beta_hm*Rel_R03
  Beta_hm_4 <- Beta_hm*Rel_R04
  

  incub <- user() ## incubation period in humans
  eip <- user() ## EIP
  inf_per <- user()  ## duration of human infectiousness
  dur_cross_prot <- user() ## duration of human cross-protection after infection
  nu <- DT/dur_cross_prot ##  waning rate of cross-protection 
  
  ## approx R0s resulting from param choices
  
  R0_1 <- kappa*kappa*Mwt*Beta_hm_mean*inf_per*Beta_mh_mean*Rel_R01/(1+delta0*eip)/delta0/R0agescale ##  eqn correct for exp distrib eip
  R0_2 <- kappa*kappa*Mwt*Beta_hm_mean*inf_per*Beta_mh_mean*Rel_R02/(1+delta0*eip)/delta0/R0agescale
  R0_3 <- kappa*kappa*Mwt*Beta_hm_mean*inf_per*Beta_mh_mean*Rel_R03/(1+delta0*eip)/delta0/R0agescale
  R0_4 <- kappa*kappa*Mwt*Beta_hm_mean*inf_per*Beta_mh_mean*Rel_R04/(1+delta0*eip)/delta0/R0agescale
  eq_FOI1 <- R0_1/init_lifespan ## in time units of years
  eq_FOI2 <- R0_2/init_lifespan
  eq_FOI3 <- R0_3/init_lifespan
  eq_FOI4 <- R0_4/init_lifespan
  
  

  ## rho params are susceptibility to infection
  ## assume no immunity beyond initial R state following each infection
  
  rho_pri <- 1
  rho_sec <- 1
  rho_tert <- 1
  rho_quart <- 1
  

  ## phi params relate to infectiousness
  
  ##  infectiousness of symptomatic cases relative to asymp
  phi_dis_enhance <- user()
  phi_ed <- phi_dis_enhance-1
  phi_scale[1:4] <- 1/(1+dis_pri[i]*phi_ed)  ##  to keep R0 definition valid
  phi_pri[1:4] <- phi_scale[i]
  phi_sec[1:4] <- phi_scale[i]
  phi_tert[1:4] <- phi_scale[i]
  phi_quart[1:4] <- phi_scale[i]
  

  ## phi_scale <- (1+0.45*phi_ed+p9p*(2+0.95*phi_ed))/(1+dis_pri*phi_ed+p9p*(2+(dis_sec+dis_tert)*phi_ed))/(1+0.45*phi_ed)  ##  to keep p9 approx accurate
  
 
  
  ################ Main model initialisation ##############
  
  ## Mosquito model first
  ## mosquito model is deterministic currently
  
  dim(init_inf_mos) <- 4
  init_inf_mos[1:4] <- user()
  
  initial(Lwt) <- Mwt*N_eq*delta0/epsilon
  initial(Mwt_S) <- Mwt*N_eq*0.9997
  initial(Mwt_E1) <- 0
  initial(Mwt_E2) <- 0
  initial(Mwt_E3) <- 0
  initial(Mwt_E4) <- 0
  initial(Mwt_I1) <- floor(Mwt*N_eq*init_inf_mos[1])
  initial(Mwt_I2) <- floor(Mwt*N_eq*init_inf_mos[2])
  initial(Mwt_I3) <- floor(Mwt*N_eq*init_inf_mos[3])
  initial(Mwt_I4) <- floor(Mwt*N_eq*init_inf_mos[4])

  Mwt_tot <- Mwt_S+Mwt_E1+Mwt_E2+Mwt_E3+Mwt_E4+Mwt_I1+Mwt_I2+Mwt_I3+Mwt_I4
  initial(R0t_1) <- 0
  initial(R0t_2) <- 0
  initial(R0t_3) <- 0
  initial(R0t_4) <- 0
  
  update(R0t_1) <- kappa*kappa*(Mwt_tot)*Beta_hm_1*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  update(R0t_2) <- kappa*kappa*(Mwt_tot)*Beta_hm_2*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  update(R0t_3) <- kappa*kappa*(Mwt_tot)*Beta_hm_3*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  update(R0t_4) <- kappa*kappa*(Mwt_tot)*Beta_hm_4*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  
  initial(out_temp) <- 0
  update(out_temp) <- temperature
  initial(out_delta) <- 0
  update(out_delta) <- delta
  initial(out_beta) <- 0
  update(out_beta) <- Beta_mh
  initial(out_Mwt) <- 0
  update(out_Mwt) <- Mwt_tot/NT
  
  Lwt_birth <- DT*gamma*Mwt_tot
  L_deathrt <- DT*sigma*((1+(Lwt/(Kc*NT))^omega))
  O_Lwt <- (1-exp(-(DT*epsilon+L_deathrt)))*(Lwt)
  Lwt_mature <- (DT*epsilon/(DT*epsilon+L_deathrt))*(O_Lwt)
  update(Lwt) <- Lwt_birth+Lwt-O_Lwt
  Mwt_FOI1 <- DT*Beta_hm_1*kappa*infectious1
  Mwt_FOI2 <- DT*Beta_hm_2*kappa*infectious2
  Mwt_FOI3 <- DT*Beta_hm_3*kappa*infectious3
  Mwt_FOI4 <- DT*Beta_hm_4*kappa*infectious4
  O_Mwt_S <- (DT*delta+Mwt_FOI1+Mwt_FOI2+Mwt_FOI3+Mwt_FOI4)*(Mwt_S)
  Mwt_inf1 <- (Mwt_FOI1/(DT*delta+Mwt_FOI1+Mwt_FOI2+Mwt_FOI3+Mwt_FOI4))*(O_Mwt_S)
  Mwt_inf2 <- (Mwt_FOI2/(DT*delta+Mwt_FOI2+Mwt_FOI3+Mwt_FOI4))*(O_Mwt_S-Mwt_inf1)
  Mwt_inf3 <- (Mwt_FOI3/(DT*delta+Mwt_FOI3+Mwt_FOI4))*(O_Mwt_S-Mwt_inf1-Mwt_inf2)
  Mwt_inf4 <- (Mwt_FOI4/(DT*delta+Mwt_FOI4))*(O_Mwt_S-Mwt_inf1-Mwt_inf2-Mwt_inf3)
  update(Mwt_S) <- Lwt_mature+Mwt_S-O_Mwt_S
  O_Mwt_E1 <- (DT*(delta+1/eip))*(Mwt_E1)
  Mwt_E1_incub <- (1/(delta*eip+1))*(O_Mwt_E1)
  update(Mwt_E1) <- Mwt_inf1+Mwt_E1-O_Mwt_E1
  O_Mwt_E2 <- (DT*(delta+1/eip))*(Mwt_E2)
  Mwt_E2_incub <- (1/(delta*eip+1))*(O_Mwt_E2)
  update(Mwt_E2) <- Mwt_inf2+Mwt_E2-O_Mwt_E2
  O_Mwt_E3 <- (DT*(delta+1/eip))*(Mwt_E3)
  Mwt_E3_incub <- (1/(delta*eip+1))*(O_Mwt_E3)
  update(Mwt_E3) <- Mwt_inf3+Mwt_E3-O_Mwt_E3
  O_Mwt_E4 <- (DT*(delta+1/eip))*(Mwt_E4)
  Mwt_E4_incub <- (1/(delta*eip+1))*(O_Mwt_E4)
  update(Mwt_E4) <- Mwt_inf4+Mwt_E4-O_Mwt_E4
  update(Mwt_I1) <- Mwt_E1_incub+Mwt_I1-(DT*delta)*(Mwt_I1)
  update(Mwt_I2) <- Mwt_E2_incub+Mwt_I2-(DT*delta)*(Mwt_I2)
  update(Mwt_I3) <- Mwt_E3_incub+Mwt_I3-(DT*delta)*(Mwt_I3)
  update(Mwt_I4) <- Mwt_E4_incub+Mwt_I4-(DT*delta)*(Mwt_I4)
  
  #### Human model ###
  
  ## state variables are vectors of dimension N_age


  initial(S[1:N_age]) <- floor(0.5+N_init_age[i]*exp(-(eq_FOI1+eq_FOI2+eq_FOI3+eq_FOI4)*mean_age[i]))

  initial(I1[1:N_age]) <- 0
  initial(I2[1:N_age]) <- 0
  initial(I3[1:N_age]) <- 0
  initial(I4[1:N_age]) <- 0
  
  initial(R1[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*exp(-(eq_FOI2+eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(R2[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*exp(-(eq_FOI1+eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(R3[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI3*mean_age[i]))*exp(-(eq_FOI1+eq_FOI2+eq_FOI4)*mean_age[i]))
  initial(R4[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI1+eq_FOI2+eq_FOI3)*mean_age[i]))

  initial(I21[1:N_age]) <- 0
  initial(I31[1:N_age]) <- 0
  initial(I41[1:N_age]) <- 0
  initial(I12[1:N_age]) <- 0
  initial(I32[1:N_age]) <- 0
  initial(I42[1:N_age]) <- 0
  initial(I13[1:N_age]) <- 0
  initial(I23[1:N_age]) <- 0
  initial(I43[1:N_age]) <- 0
  initial(I14[1:N_age]) <- 0
  initial(I24[1:N_age]) <- 0
  initial(I34[1:N_age]) <- 0
  
  initial(R12[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*exp(-(eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(R13[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*exp(-(eq_FOI2+eq_FOI4)*mean_age[i]))
  initial(R14[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI2+eq_FOI3)*mean_age[i]))
  initial(R23[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*exp(-(eq_FOI1+eq_FOI4)*mean_age[i]))
  initial(R24[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI1+eq_FOI3)*mean_age[i]))
  initial(R34[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI1+eq_FOI2)*mean_age[i]))
  
  initial(I231[1:N_age]) <- 0
  initial(I241[1:N_age]) <- 0
  initial(I341[1:N_age]) <- 0
  initial(I132[1:N_age]) <- 0
  initial(I142[1:N_age]) <- 0
  initial(I342[1:N_age]) <- 0
  initial(I123[1:N_age]) <- 0
  initial(I143[1:N_age]) <- 0
  initial(I243[1:N_age]) <- 0
  initial(I124[1:N_age]) <- 0
  initial(I134[1:N_age]) <- 0
  initial(I234[1:N_age]) <- 0
  
  initial(R123[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*exp(-eq_FOI4*mean_age[i]))
  initial(R124[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-eq_FOI3*mean_age[i]))
  initial(R134[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-eq_FOI2*mean_age[i]))
  initial(R234[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-eq_FOI1*mean_age[i]))

  
  initial(I2341[1:N_age]) <- 0
  initial(I1342[1:N_age]) <- 0
  initial(I1243[1:N_age]) <- 0
  initial(I1234[1:N_age]) <- 0
  
  initial(R1234[1:N_age]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i])))
  
  Ntotal[1:N_age] <- S[i]+I1[i]+I2[i]+I3[i]+I4[i]+R1[i]+R2[i]+R3[i]+R4[i]+I21[i]+I31[i]+I41[i]+I12[i]+I32[i]+I42[i]+I13[i]+I23[i]+I43[i]+I14[i]+I24[i]+I34[i]+R12[i]+R13[i]+R14[i]+R23[i]+R24[i]+R34[i]+I231[i]+I241[i]+I341[i]+I132[i]+I142[i]+I342[i]+I123[i]+I143[i]+I243[i]+I124[i]+I134[i]+I234[i]+R123[i]+R124[i]+R134[i]+R234[i]+I2341[i]+I1342[i]+I1243[i]+I1234[i]+R1234[i]
  NT <- sum(Ntotal)

  ## these are just derived variables
  
  initial(Ntotal_out[1:N_age]) <- N_init_age[i]
  update(Ntotal_out[1:N_age]) <- Ntotal[i]
  initial(NT_out) <- sum(N_init_age)
  update(NT_out) <- NT

  
  ## infectiousness weighted aggregated infectious compartments
  
  
  Y1[1:N_age] <- phi_pri[1]*(1+dis_pri[1]*phi_ed)*inf_1[i]+phi_sec[1]*(1+dis_sec[1]*phi_ed)*(inf_21[i]+inf_31[i]+inf_41[i])+phi_tert[1]*(1+dis_tert[1]*phi_ed)*(inf_231[i]+inf_241[i]+inf_341[i])+phi_quart[1]*(1+dis_quart[1]*phi_ed)*inf_2341[i]
  Y2[1:N_age] <- phi_pri[2]*(1+dis_pri[2]*phi_ed)*inf_2[i]+phi_sec[2]*(1+dis_sec[2]*phi_ed)*(inf_12[i]+inf_32[i]+inf_42[i])+phi_tert[2]*(1+dis_tert[2]*phi_ed)*(inf_132[i]+inf_142[i]+inf_342[i])+phi_quart[2]*(1+dis_quart[2]*phi_ed)*inf_1342[i]
  Y3[1:N_age] <- phi_pri[3]*(1+dis_pri[3]*phi_ed)*inf_3[i]+phi_sec[3]*(1+dis_sec[3]*phi_ed)*(inf_13[i]+inf_23[i]+inf_43[i])+phi_tert[3]*(1+dis_tert[3]*phi_ed)*(inf_123[i]+inf_143[i]+inf_243[i])+phi_quart[3]*(1+dis_quart[3]*phi_ed)*inf_1243[i]
  Y4[1:N_age] <- phi_pri[4]*(1+dis_pri[4]*phi_ed)*inf_4[i]+phi_sec[4]*(1+dis_sec[4]*phi_ed)*(inf_14[i]+inf_24[i]+inf_34[i])+phi_tert[4]*(1+dis_tert[4]*phi_ed)*(inf_124[i]+inf_134[i]+inf_234[i])+phi_quart[4]*(1+dis_quart[4]*phi_ed)*inf_1234[i]
  
  
  Y1T <- sum(Y1)/NT
  Y2T <- sum(Y2)/NT
  Y3T <- sum(Y3)/NT
  Y4T <- sum(Y4)/NT
  
  initial(Y1T_out) <- 0
  update(Y1T_out) <- Y1T
  initial(Y2T_out) <- 0
  update(Y2T_out) <- Y2T
  initial(Y3T_out) <- 0
  update(Y3T_out) <- Y3T
  initial(Y4T_out) <- 0
  update(Y4T_out) <- Y4T
  
  ## background force of infection (all serotypes) to keep infection going
  
  extInfRand <- extInf ## random(0,2*extInf)
  
  ## model human incubation and infectious period with aggregate serotype specific overlapping compartments
  ## rather than for every single infected state
  
  initial(exposed1) <- 0
  initial(exposed2) <- 0
  initial(exposed3) <- 0
  initial(exposed4) <- 0
  update(exposed1) <- exposed1 + Y1T - exposed1*DT/incub
  update(exposed2) <- exposed2 + Y2T - exposed2*DT/incub
  update(exposed3) <- exposed3 + Y3T - exposed3*DT/incub
  update(exposed4) <- exposed4 + Y4T - exposed4*DT/incub
  
  initial(infectious1) <- 0
  initial(infectious2) <- 0
  initial(infectious3) <- 0
  initial(infectious4) <- 0
  update(infectious1) <- infectious1*(1-DT/inf_per) + exposed1*DT/incub
  update(infectious2) <- infectious2*(1-DT/inf_per) + exposed2*DT/incub
  update(infectious3) <- infectious3*(1-DT/inf_per) + exposed3*DT/incub
  update(infectious4) <- infectious4*(1-DT/inf_per) + exposed4*DT/incub
  
  ## FOI on people (multiplied by DT)

  FOI1 <- DT*Beta_mh*kappa*Mwt_I1/NT+extInfRand
  FOI2 <- DT*Beta_mh*kappa*Mwt_I2/NT+extInfRand
  FOI3 <- DT*Beta_mh*kappa*Mwt_I3/NT+extInfRand
  FOI4 <- DT*Beta_mh*kappa*Mwt_I4/NT+extInfRand
  
  ## allow for age-related exposure
  
  FOI1a[1:N_age] <- FOIas[i]*FOI1
  FOI2a[1:N_age] <- FOIas[i]*FOI2
  FOI3a[1:N_age] <- FOIas[i]*FOI3
  FOI4a[1:N_age] <- FOIas[i]*FOI4
  
  #### derived variables to output
  
  ## output cumulative FOI
  
  initial(inc_FOI1) <- 0
  initial(inc_FOI2) <- 0
  initial(inc_FOI3) <- 0
  initial(inc_FOI4) <- 0
  update(inc_FOI1) <- inc_FOI1+FOI1
  update(inc_FOI2) <- inc_FOI2+FOI2
  update(inc_FOI3) <- inc_FOI3+FOI3
  update(inc_FOI4) <- inc_FOI4+FOI4
  initial(tot_inc_FOI) <- 0
  update(tot_inc_FOI) <- (if(NUM_YEAR_ACCUM*floor(YEAR/NUM_YEAR_ACCUM)==YEAR) 0 else tot_inc_FOI)+FOI1+FOI2+FOI3+FOI4
  initial(cur_tot_FOI) <-0
  update(cur_tot_FOI) <- FOI1+FOI2+FOI3+FOI4

  ## track and output cumulative annual human infection incidence (primary, secondary, tert+quart)
  
  infection_pri[1:N_age] <- inf_1[i]+inf_2[i]+inf_3[i]+inf_4[i]
  infection_sec[1:N_age] <- inf_21[i]+inf_31[i]+inf_41[i]+inf_12[i]+inf_32[i]+inf_42[i]+inf_13[i]+inf_23[i]+inf_43[i]+inf_14[i]+inf_24[i]+inf_34[i]
  infection_tq[1:N_age] <- inf_231[i]+inf_241[i]+inf_341[i]+inf_132[i]+inf_142[i]+inf_342[i]+inf_123[i]+inf_143[i]+inf_243[i]+inf_124[i]+inf_134[i]+inf_234[i]+inf_2341[i]+inf_1342[i]+inf_1243[i]+inf_1234[i]
  
  
  ###  outputs:weekly case incidence by age
  
  ## Not as efficient as just storing current week, but will allow entire run to be completed without interrupting
  ## incremental output would be needed for PMCMC, while normal MCMC (perhaps using summary stats) can use out_disease_array below
  

  dim(disease_sero) <- c(N_age, 4)
  disease_sero[1:N_age,1] <- dis_pri[1]*inf_1[i]+dis_sec[1]*(inf_21[i]+inf_31[i]+inf_41[i])+dis_tert[1]*(inf_231[i]+inf_241[i]+inf_341[i])+dis_quart[1]*inf_2341[i]
  disease_sero[1:N_age,2] <- dis_pri[2]*inf_2[i]+dis_sec[2]*(inf_12[i]+inf_32[i]+inf_42[i])+dis_tert[2]*(inf_132[i]+inf_142[i]+inf_342[i])+dis_quart[2]*inf_1342[i]
  disease_sero[1:N_age,3] <- dis_pri[3]*inf_3[i]+dis_sec[3]*(inf_13[i]+inf_23[i]+inf_43[i])+dis_tert[3]*(inf_123[i]+inf_143[i]+inf_243[i])+dis_quart[3]*inf_1243[i]
  disease_sero[1:N_age,4] <- dis_pri[4]*inf_4[i]+dis_sec[4]*(inf_14[i]+inf_24[i]+inf_34[i])+dis_tert[4]*(inf_124[i]+inf_134[i]+inf_234[i])+dis_quart[4]*inf_1234[i]
  dim(disease) <- N_age
  disease[1:N_age] <- sum(disease_sero[i,])
  ## dim(dis_inc) <- N_age
  ## dis_inc[1:N_age] <- disease[i]/Ntotal[i]*100000 ## incidence per 100k by age group
  disease_tot <- sum(disease[]) ## /NT*100000 ## incidence per 100k in whole population
  
  ## current week's cases (resets each week)
  dim(out_disease_cur_week) <- N_age
  initial(out_disease_cur_week[1:N_age]) <- 0.0
  update(out_disease_cur_week[1:N_age]) <- (if(cur_data_day%%7==0) 0 else out_disease_cur_week[i])+disease[i]

  # ## array of weekly cases by age (for DATA_NUM_WEEKS)
  dim(out_update_switch) <- DATA_NUM_WEEKS
  out_update_switch[1:DATA_NUM_WEEKS] <- if((cur_data_week==i) && (DATA_PERIOD_ACTIVE==1) && (cur_data_week <= DATA_NUM_WEEKS)) 1 else 0
  # dim(out_disease) <- c(DATA_NUM_WEEKS,N_age)
  # initial(out_disease[1:DATA_NUM_WEEKS,1:N_age]) <- 0.0
  # update(out_disease[1:DATA_NUM_WEEKS,1:N_age]) <- out_disease[i,j] + out_update_switch[i]*disease[j]
  # 
  # ## array of total weekly cases (for DATA_NUM_WEEKS)
  dim(out_disease_tot) <- DATA_NUM_WEEKS
  initial(out_disease_tot[1:DATA_NUM_WEEKS]) <- 0.0
  update(out_disease_tot[1:DATA_NUM_WEEKS]) <- out_disease_tot[i] + out_update_switch[i]*disease_tot
  
  ## array of total cases (over data period) by age
  dim(out_disease_age) <- N_age
  initial(out_disease_age[1:N_age]) <- 0.0
  update(out_disease_age[1:N_age]) <- out_disease_age[i] + (if((DATA_PERIOD_ACTIVE==1) && (cur_data_week <= DATA_NUM_WEEKS)) disease[i] else 0)
  
  # ## cumulative cases by ISO week and age
  dim(out_iso_update_switch) <- 53
  out_iso_update_switch[1:53] <- if((CUR_ISO_WEEK==i) && (DATA_PERIOD_ACTIVE==1) && (cur_data_week <= DATA_NUM_WEEKS)) 1 else 0
  # dim(out_disease_isowk_age) <- c(53,N_age)
  # initial(out_disease_isowk_age[1:53,1:N_age]) <- 0.0
  # update(out_disease_isowk_age[1:53,1:N_age]) <- out_disease_isowk_age[i,j] + out_iso_update_switch[i]*disease[j]
  
  ## cumulative cases by ISO week
  dim(out_disease_isowk) <- 53
  initial(out_disease_isowk[1:53]) <- 0.0
  update(out_disease_isowk[1:53]) <- out_disease_isowk[i] + out_iso_update_switch[i]*disease_tot
  
  ## cumulative cases by year
  dim(out_cases_yr_switch) <- DATA_NUM_YEARS
  out_cases_yr_switch[1:DATA_NUM_YEARS] <- if((DATA_YEAR==DATA_START_YEAR+i-1) && (DATA_PERIOD_ACTIVE==1) && (cur_data_week <= DATA_NUM_WEEKS)) 1 else 0
  dim(out_disease_tot_yr) <- DATA_NUM_YEARS
  initial(out_disease_tot_yr[1:DATA_NUM_YEARS]) <- 0
  update(out_disease_tot_yr[1:DATA_NUM_YEARS]) <- out_disease_tot_yr[i] + out_cases_yr_switch[i]*disease_tot
  
  ## seroprevalence averaged across DATA_SERO_YEAR
  dim(out_seronegative) <- N_age
  dim(seroneg) <- N_age
  seroneg[1:N_age] <- S[i]/Ntotal[i]*DT/CUR_YEAR_LENGTH
  initial(out_seronegative[1:N_age]) <- 0.0
  update(out_seronegative[1:N_age]) <- out_seronegative[i] + (if(ORDINAL_YEAR==DATA_SERO_YEAR) seroneg[i] else 0)

  ###
  ### update equations for human state variables
  ###
  
  
  O_S[1:N_age] <- rbinom((S[i]), (rho_pri*FOI1a[i]+rho_pri*FOI2a[i]+rho_pri*FOI3a[i]+rho_pri*FOI4a[i]+deathrt[i]))
  inf_1[1:N_age] <- rbinom((O_S[i]), (rho_pri*FOI1a[i]/(rho_pri*FOI1a[i]+rho_pri*FOI2a[i]+rho_pri*FOI3a[i]+rho_pri*FOI4a[i]+deathrt[i])))
  inf_2[1:N_age] <- rbinom((O_S[i]-inf_1[i]), (rho_pri*FOI2a[i]/(rho_pri*FOI2a[i]+rho_pri*FOI3a[i]+rho_pri*FOI4a[i]+deathrt[i])))
  inf_3[1:N_age] <- rbinom((O_S[i]-inf_1[i]-inf_2[i]), (rho_pri*FOI3a[i]/(rho_pri*FOI3a[i]+rho_pri*FOI4a[i]+deathrt[i])))
  inf_4[1:N_age] <- rbinom((O_S[i]-inf_1[i]-inf_2[i]-inf_3[i]), (rho_pri*FOI4a[i]/(rho_pri*FOI4a[i]+deathrt[i])))
  age_S[2:N_age_p1] <- rbinom(S[i-1]-O_S[i-1], agert[i-1])
  age_S[1] <- rbinom(10000000, births/10000000) # rpois(births)
  update(S[1:N_age]) <-  age_S[i]+S[i]-O_S[i]-age_S[i+1]
         
   O_I1[1:N_age] <- rbinom((I1[i]), (nu+deathrt[i]))
   recov_1[1:N_age] <- rbinom((O_I1[i]), (nu/(nu+deathrt[i])))
   age_I1[2:N_age_p1] <- rbinom((I1[i-1]-O_I1[i-1]+inf_1[i-1]), agert[i-1])
   age_I1[1] <- 0
   update(I1[1:N_age]) <- age_I1[i]+I1[i]-O_I1[i]+inf_1[i]-age_I1[i+1]
   
   O_I2[1:N_age] <- rbinom((I2[i]), (nu+deathrt[i]))
   recov_2[1:N_age] <- rbinom((O_I2[i]), (nu/(nu+deathrt[i])))
   age_I2[2:N_age_p1] <- rbinom((I2[i-1]-O_I2[i-1]+inf_2[i-1]), agert[i-1])
   age_I2[1] <- 0
   update(I2[1:N_age]) <- age_I2[i]+I2[i]-O_I2[i]+inf_2[i]-age_I2[i+1]
   
   O_I3[1:N_age] <- rbinom((I3[i]), (nu+deathrt[i]))
   recov_3[1:N_age] <- rbinom((O_I3[i]), (nu/(nu+deathrt[i])))
   age_I3[2:N_age_p1] <- rbinom((I3[i-1]-O_I3[i-1]+inf_3[i-1]), agert[i-1])
   age_I3[1] <- 0
   update(I3[1:N_age]) <- age_I3[i]+I3[i]-O_I3[i]+inf_3[i]-age_I3[i+1]
   
   O_I4[1:N_age] <- rbinom((I4[i]), (nu+deathrt[i]))
   recov_4[1:N_age] <- rbinom((O_I4[i]), (nu/(nu+deathrt[i])))
   age_I4[2:N_age_p1] <- rbinom((I4[i-1]-O_I4[i-1]+inf_4[i-1]), agert[i-1])
   age_I4[1] <- 0
   update(I4[1:N_age]) <- age_I4[i]+I4[i]-O_I4[i]+inf_4[i]-age_I4[i+1]
   
   O_R1[1:N_age] <- rbinom((R1[i]), (rho_sec*FOI2a[i]+rho_sec*FOI3a[i]+rho_sec*FOI4a[i]+deathrt[i]))
   inf_12[1:N_age] <- rbinom((O_R1[i]), (rho_sec*FOI2a[i]/(rho_sec*FOI2a[i]+rho_sec*FOI3a[i]+rho_sec*FOI4a[i]+deathrt[i])))
   inf_13[1:N_age] <- rbinom((O_R1[i]-inf_12[i]), (rho_sec*FOI3a[i]/(rho_sec*FOI3a[i]+rho_sec*FOI4a[i]+deathrt[i])))
   inf_14[1:N_age] <- rbinom((O_R1[i]-inf_12[i]-inf_13[i]), (rho_sec*FOI4a[i]/(rho_sec*FOI4a[i]+deathrt[i])))
   age_R1[2:N_age_p1] <- rbinom((R1[i-1]-O_R1[i-1]+recov_1[i-1]), agert[i-1])
   age_R1[1] <- 0
   update(R1[1:N_age]) <- age_R1[i]+R1[i]-O_R1[i]+recov_1[i]-age_R1[i+1]
   
   O_R2[1:N_age] <- rbinom((R2[i]), (rho_sec*FOI1a[i]+rho_sec*FOI3a[i]+rho_sec*FOI4a[i]+deathrt[i]))
   inf_21[1:N_age] <- rbinom((O_R2[i]), (rho_sec*FOI1a[i]/(rho_sec*FOI1a[i]+rho_sec*FOI3a[i]+rho_sec*FOI4a[i]+deathrt[i])))
   inf_23[1:N_age] <- rbinom((O_R2[i]-inf_21[i]), (rho_sec*FOI3a[i]/(rho_sec*FOI3a[i]+rho_sec*FOI4a[i]+deathrt[i])))
   inf_24[1:N_age] <- rbinom((O_R2[i]-inf_21[i]-inf_23[i]), (rho_sec*FOI4a[i]/(rho_sec*FOI4a[i]+deathrt[i])))
   age_R2[2:N_age_p1] <- rbinom((R2[i-1]-O_R2[i-1]+recov_2[i-1]), agert[i-1])
   age_R2[1] <- 0 
   update(R2[1:N_age]) <- age_R2[i]+R2[i]-O_R2[i]+recov_2[i]-age_R2[i+1]
   
   O_R3[1:N_age] <- rbinom((R3[i]), (rho_sec*FOI1a[i]+rho_sec*FOI2a[i]+rho_sec*FOI4a[i]+deathrt[i]))
   inf_31[1:N_age] <- rbinom((O_R3[i]), (rho_sec*FOI1a[i]/(rho_sec*FOI1a[i]+rho_sec*FOI2a[i]+rho_sec*FOI4a[i]+deathrt[i])))
   inf_32[1:N_age] <- rbinom((O_R3[i]-inf_31[i]), (rho_sec*FOI2a[i]/(rho_sec*FOI2a[i]+rho_sec*FOI4a[i]+deathrt[i])))
   inf_34[1:N_age] <- rbinom((O_R3[i]-inf_31[i]-inf_32[i]), (rho_sec*FOI4a[i]/(rho_sec*FOI4a[i]+deathrt[i])))
   age_R3[2:N_age_p1] <- rbinom((R3[i-1]-O_R3[i-1]+recov_3[i-1]), agert[i-1])
   age_R3[1] <- 0 
   update(R3[1:N_age]) <- age_R3[i]+R3[i]-O_R3[i]+recov_3[i]-age_R3[i+1]
   
   O_R4[1:N_age] <- rbinom((R4[i]), (rho_sec*FOI1a[i]+rho_sec*FOI2a[i]+rho_sec*FOI3a[i]+deathrt[i]))
   inf_41[1:N_age] <- rbinom((O_R4[i]), (rho_sec*FOI1a[i]/(rho_sec*FOI1a[i]+rho_sec*FOI2a[i]+rho_sec*FOI3a[i]+deathrt[i])))
   inf_42[1:N_age] <- rbinom((O_R4[i]-inf_41[i]), (rho_sec*FOI2a[i]/(rho_sec*FOI2a[i]+rho_sec*FOI3a[i]+deathrt[i])))
   inf_43[1:N_age] <- rbinom((O_R4[i]-inf_41[i]-inf_42[i]), (rho_sec*FOI3a[i]/(rho_sec*FOI3a[i]+deathrt[i])))
   age_R4[2:N_age_p1] <- rbinom((R4[i-1]-O_R4[i-1]+recov_4[i-1]), agert[i-1])
   age_R4[1] <- 0
   update(R4[1:N_age]) <- age_R4[i]+R4[i]-O_R4[i]+recov_4[i]-age_R4[i+1]
   
   O_I12[1:N_age] <- rbinom((I12[i]), (nu+deathrt[i]))
   recov_12[1:N_age] <- rbinom((O_I12[i]), (nu/(nu+deathrt[i])))
   age_I12[2:N_age_p1] <- rbinom((I12[i-1]-O_I12[i-1]+inf_12[i-1]), agert[i-1])
   age_I12[1] <- 0
   update(I12[1:N_age]) <- age_I12[i]+I12[i]-O_I12[i]+inf_12[i]-age_I12[i+1]
   
   O_I13[1:N_age] <- rbinom((I13[i]), (nu+deathrt[i]))
   recov_13[1:N_age] <- rbinom((O_I13[i]), (nu/(nu+deathrt[i])))
   age_I13[2:N_age_p1] <- rbinom((I13[i-1]-O_I13[i-1]+inf_13[i-1]), agert[i-1])
   age_I13[1] <- 0
   update(I13[1:N_age]) <- age_I13[i]+I13[i]-O_I13[i]+inf_13[i]-age_I13[i+1]
   
   O_I14[1:N_age] <- rbinom((I14[i]), (nu+deathrt[i]))
   recov_14[1:N_age] <- rbinom((O_I14[i]), (nu/(nu+deathrt[i])))
   age_I14[2:N_age_p1] <- rbinom((I14[i-1]-O_I14[i-1]+inf_14[i-1]), agert[i-1])
   age_I14[1] <- 0
   update(I14[1:N_age]) <- age_I14[i]+I14[i]-O_I14[i]+inf_14[i]-age_I14[i+1]
   
   O_I21[1:N_age] <- rbinom((I21[i]), (nu+deathrt[i]))
   recov_21[1:N_age] <- rbinom((O_I21[i]), (nu/(nu+deathrt[i])))
   age_I21[2:N_age_p1] <- rbinom((I21[i-1]-O_I21[i-1]+inf_21[i-1]), agert[i-1])
   age_I21[1] <- 0
   update(I21[1:N_age]) <- age_I21[i]+I21[i]-O_I21[i]+inf_21[i]-age_I21[i+1]
   
   O_I23[1:N_age] <- rbinom((I23[i]), (nu+deathrt[i]))
   recov_23[1:N_age] <- rbinom((O_I23[i]), (nu/(nu+deathrt[i])))
   age_I23[2:N_age_p1] <- rbinom((I23[i-1]-O_I23[i-1]+inf_23[i-1]), agert[i-1])
   age_I23[1] <- 0
   update(I23[1:N_age]) <- age_I23[i]+I23[i]-O_I23[i]+inf_23[i]-age_I23[i+1]
   
   O_I24[1:N_age] <- rbinom((I24[i]), (nu+deathrt[i]))
   recov_24[1:N_age] <- rbinom((O_I24[i]), (nu/(nu+deathrt[i])))
   age_I24[2:N_age_p1] <- rbinom((I24[i-1]-O_I24[i-1]+inf_24[i-1]), agert[i-1])
   age_I24[1] <- 0
   update(I24[1:N_age]) <- age_I24[i]+I24[i]-O_I24[i]+inf_24[i]-age_I24[i+1]
   
   O_I31[1:N_age] <- rbinom((I31[i]), (nu+deathrt[i]))
   recov_31[1:N_age] <- rbinom((O_I31[i]), (nu/(nu+deathrt[i])))
   age_I31[2:N_age_p1] <- rbinom((I31[i-1]-O_I31[i-1]+inf_31[i-1]), agert[i-1])
   age_I31[1] <- 0
   update(I31[1:N_age]) <- age_I31[i]+I31[i]-O_I31[i]+inf_31[i]-age_I31[i+1]
   
   O_I32[1:N_age] <- rbinom((I32[i]), (nu+deathrt[i]))
   recov_32[1:N_age] <- rbinom((O_I32[i]), (nu/(nu+deathrt[i])))
   age_I32[2:N_age_p1] <- rbinom((I32[i-1]-O_I32[i-1]+inf_32[i-1]), agert[i-1])
   age_I32[1] <- 0
   update(I32[1:N_age]) <- age_I32[i]+I32[i]-O_I32[i]+inf_32[i]-age_I32[i+1]
   
   O_I34[1:N_age] <- rbinom((I34[i]), (nu+deathrt[i]))
   recov_34[1:N_age] <- rbinom((O_I34[i]), (nu/(nu+deathrt[i])))
   age_I34[2:N_age_p1] <- rbinom((I34[i-1]-O_I34[i-1]+inf_34[i-1]), agert[i-1])
   age_I34[1] <- 0
   update(I34[1:N_age]) <- age_I34[i]+I34[i]-O_I34[i]+inf_34[i]-age_I34[i+1]
   
   O_I41[1:N_age] <- rbinom((I41[i]), (nu+deathrt[i]))
   recov_41[1:N_age] <- rbinom((O_I41[i]), (nu/(nu+deathrt[i])))
   age_I41[2:N_age_p1] <- rbinom((I41[i-1]-O_I41[i-1]+inf_41[i-1]), agert[i-1])
   age_I41[1] <- 0 
   update(I41[1:N_age]) <- age_I41[i]+I41[i]-O_I41[i]+inf_41[i]-age_I41[i+1]
   
   O_I42[1:N_age] <- rbinom((I42[i]), (nu+deathrt[i]))
   recov_42[1:N_age] <- rbinom((O_I42[i]), (nu/(nu+deathrt[i])))
   age_I42[2:N_age_p1] <- rbinom((I42[i-1]-O_I42[i-1]+inf_42[i-1]), agert[i-1])
   age_I42[1] <- 0
   update(I42[1:N_age]) <- age_I42[i]+I42[i]-O_I42[i]+inf_42[i]-age_I42[i+1]
   
   O_I43[1:N_age] <- rbinom((I43[i]), (nu+deathrt[i]))
   recov_43[1:N_age] <- rbinom((O_I43[i]), (nu/(nu+deathrt[i])))
   age_I43[2:N_age_p1] <- rbinom((I43[i-1]-O_I43[i-1]+inf_43[i-1]), agert[i-1])
   age_I43[1] <- 0
   update(I43[1:N_age]) <- age_I43[i]+I43[i]-O_I43[i]+inf_43[i]-age_I43[i+1]
   
   O_R12[1:N_age] <- rbinom((R12[i]), (rho_tert*FOI3a[i]+rho_tert*FOI4a[i]+deathrt[i]))
   inf_123[1:N_age] <- rbinom((O_R12[i]), (rho_tert*FOI3a[i]/(rho_tert*FOI3a[i]+rho_tert*FOI4a[i]+deathrt[i])))
   inf_124[1:N_age] <- rbinom((O_R12[i]-inf_123[i]), (rho_tert*FOI4a[i]/(rho_tert*FOI4a[i]+deathrt[i])))
   age_R12[2:N_age_p1] <- rbinom((R12[i-1]-O_R12[i-1]+recov_12[i-1]+recov_21[i-1]), agert[i-1])
   age_R12[1] <- 0
   update(R12[1:N_age]) <- age_R12[i]+R12[i]-O_R12[i]+recov_12[i]+recov_21[i]-age_R12[i+1]
   
   O_R13[1:N_age] <- rbinom((R13[i]), (rho_tert*FOI2a[i]+rho_tert*FOI4a[i]+deathrt[i]))
   inf_132[1:N_age] <- rbinom((O_R13[i]), (rho_tert*FOI2a[i]/(rho_tert*FOI2a[i]+rho_tert*FOI4a[i]+deathrt[i])))
   inf_134[1:N_age] <- rbinom((O_R13[i]-inf_132[i]), (rho_tert*FOI4a[i]/(rho_tert*FOI4a[i]+deathrt[i])))
   age_R13[2:N_age_p1] <- rbinom((R13[i-1]-O_R13[i-1]+recov_13[i-1]+recov_31[i-1]), agert[i-1])
   age_R13[1] <- 0
   update(R13[1:N_age]) <- age_R13[i]+R13[i]-O_R13[i]+recov_13[i]+recov_31[i]-age_R13[i+1]
   
   O_R14[1:N_age] <- rbinom((R14[i]), (rho_tert*FOI2a[i]+rho_tert*FOI3a[i]+deathrt[i]))
   inf_142[1:N_age] <- rbinom((O_R14[i]), (rho_tert*FOI2a[i]/(rho_tert*FOI2a[i]+rho_tert*FOI3a[i]+deathrt[i])))
   inf_143[1:N_age] <- rbinom((O_R14[i]-inf_142[i]), (rho_tert*FOI3a[i]/(rho_tert*FOI3a[i]+deathrt[i])))
   age_R14[2:N_age_p1] <- rbinom((R14[i-1]-O_R14[i-1]+recov_14[i-1]+recov_41[i-1]), agert[i-1])
   age_R14[1] <- 0
   update(R14[1:N_age]) <- age_R14[i]+R14[i]-O_R14[i]+recov_14[i]+recov_41[i]-age_R14[i+1]
   
   O_R23[1:N_age] <- rbinom((R23[i]), (rho_tert*FOI1a[i]+rho_tert*FOI4a[i]+deathrt[i]))
   inf_231[1:N_age] <- rbinom((O_R23[i]), (rho_tert*FOI1a[i]/(rho_tert*FOI1a[i]+rho_tert*FOI4a[i]+deathrt[i])))
   inf_234[1:N_age] <- rbinom((O_R23[i]-inf_231[i]), (rho_tert*FOI4a[i]/(rho_tert*FOI4a[i]+deathrt[i])))
   age_R23[2:N_age_p1] <- rbinom((R23[i-1]-O_R23[i-1]+recov_23[i-1]+recov_32[i-1]), agert[i-1])
   age_R23[1] <- 0
   update(R23[1:N_age]) <- age_R23[i]+R23[i]-O_R23[i]+recov_23[i]+recov_32[i]-age_R23[i+1]
   
   O_R24[1:N_age] <- rbinom((R24[i]), (rho_tert*FOI1a[i]+rho_tert*FOI3a[i]+deathrt[i]))
   inf_241[1:N_age] <- rbinom((O_R24[i]), (rho_tert*FOI1a[i]/(rho_tert*FOI1a[i]+rho_tert*FOI3a[i]+deathrt[i])))
   inf_243[1:N_age] <- rbinom((O_R24[i]-inf_241[i]), (rho_tert*FOI3a[i]/(rho_tert*FOI3a[i]+deathrt[i])))
   age_R24[2:N_age_p1] <- rbinom((R24[i-1]-O_R24[i-1]+recov_24[i-1]+recov_42[i-1]), agert[i-1])
   age_R24[1] <- 0
   update(R24[1:N_age]) <- age_R24[i]+R24[i]-O_R24[i]+recov_24[i]+recov_42[i]-age_R24[i+1]
   
   O_R34[1:N_age] <- rbinom((R34[i]), (rho_tert*FOI1a[i]+rho_tert*FOI2a[i]+deathrt[i]))
   inf_341[1:N_age] <- rbinom((O_R34[i]), (rho_tert*FOI1a[i]/(rho_tert*FOI1a[i]+rho_tert*FOI2a[i]+deathrt[i])))
   inf_342[1:N_age] <- rbinom((O_R34[i]-inf_341[i]), (rho_tert*FOI2a[i]/(rho_tert*FOI2a[i]+deathrt[i])))
   age_R34[2:N_age_p1] <- rbinom((R34[i-1]-O_R34[i-1]+recov_34[i-1]+recov_43[i-1]), agert[i-1])
   age_R34[1] <- 0
   update(R34[1:N_age]) <- age_R34[i]+R34[i]-O_R34[i]+recov_34[i]+recov_43[i]-age_R34[i+1]
   
   O_I123[1:N_age] <- rbinom((I123[i]), (nu+deathrt[i]))
   recov_123[1:N_age] <- rbinom((O_I123[i]), (nu/(nu+deathrt[i])))
   age_I123[2:N_age_p1] <- rbinom((I123[i-1]-O_I123[i-1]+inf_123[i-1]), agert[i-1])
   age_I123[1] <- 0
   update(I123[1:N_age]) <- age_I123[i]+I123[i]-O_I123[i]+inf_123[i]-age_I123[i+1]
   
   O_I124[1:N_age] <- rbinom((I124[i]), (nu+deathrt[i]))
   recov_124[1:N_age] <- rbinom((O_I124[i]), (nu/(nu+deathrt[i])))
   age_I124[2:N_age_p1] <- rbinom((I124[i-1]-O_I124[i-1]+inf_124[i-1]), agert[i-1])
   age_I124[1] <- 0 
   update(I124[1:N_age]) <- age_I124[i]+I124[i]-O_I124[i]+inf_124[i]-age_I124[i+1]
   
   O_I132[1:N_age] <- rbinom((I132[i]), (nu+deathrt[i]))
   recov_132[1:N_age] <- rbinom((O_I132[i]), (nu/(nu+deathrt[i])))
   age_I132[2:N_age_p1] <- rbinom((I132[i-1]-O_I132[i-1]+inf_132[i-1]), agert[i-1])
   age_I132[1] <- 0
   update(I132[1:N_age]) <- age_I132[i]+I132[i]-O_I132[i]+inf_132[i]-age_I132[i+1]
   
   O_I134[1:N_age] <- rbinom((I134[i]), (nu+deathrt[i]))
   recov_134[1:N_age] <- rbinom((O_I134[i]), (nu/(nu+deathrt[i])))
   age_I134[2:N_age_p1] <- rbinom((I134[i-1]-O_I134[i-1]+inf_134[i-1]), agert[i-1])
   age_I134[1] <- 0 
   update(I134[1:N_age]) <- age_I134[i]+I134[i]-O_I134[i]+inf_134[i]-age_I134[i+1]
   
   O_I142[1:N_age] <- rbinom((I142[i]), (nu+deathrt[i]))
   recov_142[1:N_age] <- rbinom((O_I142[i]), (nu/(nu+deathrt[i])))
   age_I142[2:N_age_p1] <- rbinom((I142[i-1]-O_I142[i-1]+inf_142[i-1]), agert[i-1])
   age_I142[1] <- 0
   update(I142[1:N_age]) <- age_I142[i]+I142[i]-O_I142[i]+inf_142[i]-age_I142[i+1]
   
   O_I143[1:N_age] <- rbinom((I143[i]), (nu+deathrt[i]))
   recov_143[1:N_age] <- rbinom((O_I143[i]), (nu/(nu+deathrt[i])))
   age_I143[2:N_age_p1] <- rbinom((I143[i-1]-O_I143[i-1]+inf_143[i-1]), agert[i-1])
   age_I143[1] <- 0 
   update(I143[1:N_age]) <- age_I143[i]+I143[i]-O_I143[i]+inf_143[i]-age_I143[i+1]
   
   O_I231[1:N_age] <- rbinom((I231[i]), (nu+deathrt[i]))
   recov_231[1:N_age] <- rbinom((O_I231[i]), (nu/(nu+deathrt[i])))
   age_I231[2:N_age_p1] <- rbinom((I231[i-1]-O_I231[i-1]+inf_231[i-1]), agert[i-1])
   age_I231[1] <- 0
   update(I231[1:N_age]) <- age_I231[i]+I231[i]-O_I231[i]+inf_231[i]-age_I231[i+1]
   
   O_I234[1:N_age] <- rbinom((I234[i]), (nu+deathrt[i]))
   recov_234[1:N_age] <- rbinom((O_I234[i]), (nu/(nu+deathrt[i])))
   age_I234[2:N_age_p1] <- rbinom((I234[i-1]-O_I234[i-1]+inf_234[i-1]), agert[i-1])
   age_I234[1] <- 0
   update(I234[1:N_age]) <- age_I234[i]+I234[i]-O_I234[i]+inf_234[i]-age_I234[i+1]
   
   O_I241[1:N_age] <- rbinom((I241[i]), (nu+deathrt[i]))
   recov_241[1:N_age] <- rbinom((O_I241[i]), (nu/(nu+deathrt[i])))
   age_I241[2:N_age_p1] <- rbinom((I241[i-1]-O_I241[i-1]+inf_241[i-1]), agert[i-1])
   age_I241[1] <- 0
   update(I241[1:N_age]) <- age_I241[i]+I241[i]-O_I241[i]+inf_241[i]-age_I241[i+1]
   
   O_I243[1:N_age] <- rbinom((I243[i]), (nu+deathrt[i]))
   recov_243[1:N_age] <- rbinom((O_I243[i]), (nu/(nu+deathrt[i])))
   age_I243[2:N_age_p1] <- rbinom((I243[i-1]-O_I243[i-1]+inf_243[i-1]), agert[i-1])
   age_I243[1] <- 0
   update(I243[1:N_age]) <- age_I243[i]+I243[i]-O_I243[i]+inf_243[i]-age_I243[i+1]
   
   O_I341[1:N_age] <- rbinom((I341[i]), (nu+deathrt[i]))
   recov_341[1:N_age] <- rbinom((O_I341[i]), (nu/(nu+deathrt[i])))
   age_I341[2:N_age_p1] <- rbinom((I341[i-1]-O_I341[i-1]+inf_341[i-1]), agert[i-1])
   age_I341[1] <- 0
   update(I341[1:N_age]) <- age_I341[i]+I341[i]-O_I341[i]+inf_341[i]-age_I341[i+1]
   
   O_I342[1:N_age] <- rbinom((I342[i]), (nu+deathrt[i]))
   recov_342[1:N_age] <- rbinom((O_I342[i]), (nu/(nu+deathrt[i])))
   age_I342[2:N_age_p1] <- rbinom((I342[i-1]-O_I342[i-1]+inf_342[i-1]), agert[i-1])
   age_I342[1] <- 0 
   update(I342[1:N_age]) <- age_I342[i]+I342[i]-O_I342[i]+inf_342[i]-age_I342[i+1]
   
   O_R123[1:N_age] <- rbinom((R123[i]), (rho_quart*FOI4a[i]+deathrt[i]))
   inf_1234[1:N_age] <- rbinom((O_R123[i]), (rho_quart*FOI4a[i]/(rho_quart*FOI4a[i]+deathrt[i])))
   age_R123[2:N_age_p1] <- rbinom((R123[i-1]-O_R123[i-1]+recov_123[i-1]+recov_132[i-1]+recov_231[i-1]), agert[i-1])
   age_R123[1] <- 0
   update(R123[1:N_age]) <- age_R123[i]+R123[i]-O_R123[i]+recov_123[i]+recov_132[i]+recov_231[i]-age_R123[i+1]
   
   O_R124[1:N_age] <- rbinom((R124[i]), (rho_quart*FOI3a[i]+deathrt[i]))
   inf_1243[1:N_age] <- rbinom((O_R124[i]), (rho_quart*FOI3a[i]/(rho_quart*FOI3a[i]+deathrt[i])))
   age_R124[2:N_age_p1] <- rbinom((R124[i-1]-O_R124[i-1]+recov_124[i-1]+recov_142[i-1]+recov_241[i-1]), agert[i-1])
   age_R124[1] <- 0
   update(R124[1:N_age]) <- age_R124[i]+R124[i]-O_R124[i]+recov_124[i]+recov_142[i]+recov_241[i]-age_R124[i+1]
   
   O_R134[1:N_age] <- rbinom((R134[i]), (rho_quart*FOI2a[i]+deathrt[i]))
   inf_1342[1:N_age] <- rbinom((O_R134[i]), (rho_quart*FOI2a[i]/(rho_quart*FOI2a[i]+deathrt[i])))
   age_R134[2:N_age_p1] <- rbinom((R134[i-1]-O_R134[i-1]+recov_134[i-1]+recov_143[i-1]+recov_341[i-1]), agert[i-1])
   age_R134[1] <- 0
   update(R134[1:N_age]) <- age_R134[i]+R134[i]-O_R134[i]+recov_134[i]+recov_143[i]+recov_341[i]-age_R134[i+1]
   
   O_R234[1:N_age] <- rbinom((R234[i]), (rho_quart*FOI1a[i]+deathrt[i]))
   inf_2341[1:N_age] <- rbinom((O_R234[i]), (rho_quart*FOI1a[i]/(rho_quart*FOI1a[i]+deathrt[i])))
   age_R234[2:N_age_p1] <- rbinom((R234[i-1]-O_R234[i-1]+recov_234[i-1]+recov_243[i-1]+recov_342[i-1]), agert[i-1])
   age_R234[1] <- 0
   update(R234[1:N_age]) <- age_R234[i]+R234[i]-O_R234[i]+recov_234[i]+recov_243[i]+recov_342[i]-age_R234[i+1]
   
   O_I1234[1:N_age] <- rbinom((I1234[i]), (nu+deathrt[i]))
   recov_1234[1:N_age] <- rbinom((O_I1234[i]), (nu/(nu+deathrt[i])))
   age_I1234[2:N_age_p1] <- rbinom((I1234[i-1]-O_I1234[i-1]+inf_1234[i-1]), agert[i-1])
   age_I1234[1] <- 0 
   update(I1234[1:N_age]) <- age_I1234[i]+I1234[i]-O_I1234[i]+inf_1234[i]-age_I1234[i+1]
   
   O_I1243[1:N_age] <- rbinom((I1243[i]), (nu+deathrt[i]))
   recov_1243[1:N_age] <- rbinom((O_I1243[i]), (nu/(nu+deathrt[i])))
   age_I1243[2:N_age_p1] <- rbinom((I1243[i-1]-O_I1243[i-1]+inf_1243[i-1]), agert[i-1])
   age_I1243[1] <- 0
   update(I1243[1:N_age]) <- age_I1243[i]+I1243[i]-O_I1243[i]+inf_1243[i]-age_I1243[i+1]
   
   O_I1342[1:N_age] <- rbinom((I1342[i]), (nu+deathrt[i]))
   recov_1342[1:N_age] <- rbinom((O_I1342[i]), (nu/(nu+deathrt[i])))
   age_I1342[2:N_age_p1] <- rbinom((I1342[i-1]-O_I1342[i-1]+inf_1342[i-1]), agert[i-1])
   age_I1342[1] <- 0
   update(I1342[1:N_age]) <- age_I1342[i]+I1342[i]-O_I1342[i]+inf_1342[i]-age_I1342[i+1]
   
   O_I2341[1:N_age] <- rbinom((I2341[i]), (nu+deathrt[i]))
   recov_2341[1:N_age] <- rbinom((O_I2341[i]), (nu/(nu+deathrt[i])))
   age_I2341[2:N_age_p1] <- rbinom((I2341[i-1]-O_I2341[i-1]+inf_2341[i-1]), agert[i-1])
   age_I2341[1] <- 0
   update(I2341[1:N_age]) <- age_I2341[i]+I2341[i]-O_I2341[i]+inf_2341[i]-age_I2341[i+1]
   
   O_R1234[1:N_age] <- rbinom((R1234[i]), (deathrt[i]))
   age_R1234[2:N_age_p1] <- rbinom((R1234[i-1]-O_R1234[i-1]+recov_1234[i-1]+recov_1243[i-1]+recov_1342[i-1]+recov_2341[i-1]), agert[i-1])
   age_R1234[1] <- 0
   update(R1234[1:N_age]) <- age_R1234[i]+R1234[i]-O_R1234[i]+recov_1234[i]+recov_1243[i]+recov_1342[i]+recov_2341[i]-age_R1234[i+1]
         
 
## array dimension

  dim(age_removal_d) <- c(151,N_age_p1)
  dim(age_rate_d) <-  c(151,N_age_p1)
  dim(pop_size_d) <-  N_age_p1
  dim(births_d) <- c(151,2)
  dim(life_expec_d) <-  c(151,N_age_p1)
  
  dim(climate_d) <- c(NTP,2) ## col 1 = temp, col 2 = rainfall, daily

  dim(dis_pri) <- 4
  dim(dis_sec) <- 4
  dim(dis_tert) <- 4
  dim(dis_quart) <- 4
  
  dim(sdis_pri) <- 4
  dim(sdis_sec) <- 4
  dim(sdis_tert) <- 4
  dim(sdis_quart) <- 4
  
  dim(agec) <- N_age
  
  dim(death) <- N_age
  dim(ageb) <- N_age_p1
  dim(mean_age) <- N_age
  dim(init_life_expec) <- N_age
  dim(deathrt) <- N_age
  dim(cur_age_rate) <- N_age
  dim(agert) <- N_age
  dim(N_init_age0) <- N_age
  dim(N_init_age) <- N_age
  dim(life_expec) <- N_age
  dim(FOIas) <- N_age
  dim(Ntotal_out) <- N_age
  dim(suscinitpop) <- N_age


  dim(phi_scale) <- 4
  dim(phi_pri) <- 4
  dim(phi_sec) <- 4
  dim(phi_tert) <- 4
  dim(phi_quart) <- 4
  
  
  dim(S) <- N_age
  
  dim(I1) <- N_age
  dim(I2) <- N_age
  dim(I3) <- N_age
  dim(I4) <- N_age
  dim(I21) <- N_age
  dim(I31) <- N_age
  dim(I41) <- N_age
  dim(I12) <- N_age
  dim(I32) <- N_age
  dim(I42) <- N_age
  dim(I13) <- N_age
  dim(I23) <- N_age
  dim(I43) <- N_age
  dim(I14) <- N_age
  dim(I24) <- N_age
  dim(I34) <- N_age
  dim(I231) <- N_age
  dim(I241) <- N_age
  dim(I341) <- N_age
  dim(I132) <- N_age
  dim(I142) <- N_age
  dim(I342) <- N_age
  dim(I123) <- N_age
  dim(I143) <- N_age
  dim(I243) <- N_age
  dim(I124) <- N_age
  dim(I134) <- N_age
  dim(I234) <- N_age
  dim(I2341) <- N_age
  dim(I1342) <- N_age
  dim(I1243) <- N_age
  dim(I1234) <- N_age
  
  dim(R1) <- N_age
  dim(R2) <- N_age
  dim(R3) <- N_age
  dim(R4) <- N_age
  dim(R12) <- N_age
  dim(R13) <- N_age
  dim(R14) <- N_age
  dim(R23) <- N_age
  dim(R24) <- N_age
  dim(R34) <- N_age
  dim(R123) <- N_age
  dim(R124) <- N_age
  dim(R134) <- N_age
  dim(R234) <- N_age
  dim(R1234) <- N_age
  
  dim(Ntotal) <- N_age

  dim(Y1) <- N_age
  dim(Y2) <- N_age
  dim(Y3) <- N_age
  dim(Y4) <- N_age
  
  dim(FOI1a) <- N_age
  dim(FOI2a) <- N_age
  dim(FOI3a) <- N_age
  dim(FOI4a) <- N_age
  
  dim(infection_pri) <- N_age
  dim(infection_sec) <- N_age
  dim(infection_tq) <- N_age

  dim(inf_1) <- N_age
  dim(inf_2) <- N_age
  dim(inf_3) <- N_age
  dim(inf_4) <- N_age
  dim(inf_21) <- N_age
  dim(inf_31) <- N_age
  dim(inf_41) <- N_age
  dim(inf_12) <- N_age
  dim(inf_32) <- N_age
  dim(inf_42) <- N_age
  dim(inf_13) <- N_age
  dim(inf_23) <- N_age
  dim(inf_43) <- N_age
  dim(inf_14) <- N_age
  dim(inf_24) <- N_age
  dim(inf_34) <- N_age
  dim(inf_231) <- N_age
  dim(inf_241) <- N_age
  dim(inf_341) <- N_age
  dim(inf_132) <- N_age
  dim(inf_142) <- N_age
  dim(inf_342) <- N_age
  dim(inf_123) <- N_age
  dim(inf_143) <- N_age
  dim(inf_243) <- N_age
  dim(inf_124) <- N_age
  dim(inf_134) <- N_age
  dim(inf_234) <- N_age
  dim(inf_2341) <- N_age
  dim(inf_1342) <- N_age
  dim(inf_1243) <- N_age
  dim(inf_1234) <- N_age
  
  dim(recov_1) <- N_age
  dim(recov_2) <- N_age
  dim(recov_3) <- N_age
  dim(recov_4) <- N_age
  dim(recov_21) <- N_age
  dim(recov_31) <- N_age
  dim(recov_41) <- N_age
  dim(recov_12) <- N_age
  dim(recov_32) <- N_age
  dim(recov_42) <- N_age
  dim(recov_13) <- N_age
  dim(recov_23) <- N_age
  dim(recov_43) <- N_age
  dim(recov_14) <- N_age
  dim(recov_24) <- N_age
  dim(recov_34) <- N_age
  dim(recov_231) <- N_age
  dim(recov_241) <- N_age
  dim(recov_341) <- N_age
  dim(recov_132) <- N_age
  dim(recov_142) <- N_age
  dim(recov_342) <- N_age
  dim(recov_123) <- N_age
  dim(recov_143) <- N_age
  dim(recov_243) <- N_age
  dim(recov_124) <- N_age
  dim(recov_134) <- N_age
  dim(recov_234) <- N_age
  dim(recov_2341) <- N_age
  dim(recov_1342) <- N_age
  dim(recov_1243) <- N_age
  dim(recov_1234) <- N_age
  
  dim(age_S) <- N_age_p1
  
  dim(age_I1) <- N_age_p1
  dim(age_I2) <- N_age_p1
  dim(age_I3) <- N_age_p1
  dim(age_I4) <- N_age_p1
  dim(age_I21) <- N_age_p1
  dim(age_I31) <- N_age_p1
  dim(age_I41) <- N_age_p1
  dim(age_I12) <- N_age_p1
  dim(age_I32) <- N_age_p1
  dim(age_I42) <- N_age_p1
  dim(age_I13) <- N_age_p1
  dim(age_I23) <- N_age_p1
  dim(age_I43) <- N_age_p1
  dim(age_I14) <- N_age_p1
  dim(age_I24) <- N_age_p1
  dim(age_I34) <- N_age_p1
  dim(age_I231) <- N_age_p1
  dim(age_I241) <- N_age_p1
  dim(age_I341) <- N_age_p1
  dim(age_I132) <- N_age_p1
  dim(age_I142) <- N_age_p1
  dim(age_I342) <- N_age_p1
  dim(age_I123) <- N_age_p1
  dim(age_I143) <- N_age_p1
  dim(age_I243) <- N_age_p1
  dim(age_I124) <- N_age_p1
  dim(age_I134) <- N_age_p1
  dim(age_I234) <- N_age_p1
  dim(age_I2341) <- N_age_p1
  dim(age_I1342) <- N_age_p1
  dim(age_I1243) <- N_age_p1
  dim(age_I1234) <- N_age_p1
  
  dim(age_R1) <- N_age_p1
  dim(age_R2) <- N_age_p1
  dim(age_R3) <- N_age_p1
  dim(age_R4) <- N_age_p1
  dim(age_R12) <- N_age_p1
  dim(age_R13) <- N_age_p1
  dim(age_R14) <- N_age_p1
  dim(age_R23) <- N_age_p1
  dim(age_R24) <- N_age_p1
  dim(age_R34) <- N_age_p1
  dim(age_R123) <- N_age_p1
  dim(age_R124) <- N_age_p1
  dim(age_R134) <- N_age_p1
  dim(age_R234) <- N_age_p1
  dim(age_R1234) <- N_age_p1
  

  dim(O_S) <- N_age
  
  dim(O_I1) <- N_age
  dim(O_I2) <- N_age
  dim(O_I3) <- N_age
  dim(O_I4) <- N_age
  dim(O_I21) <- N_age
  dim(O_I31) <- N_age
  dim(O_I41) <- N_age
  dim(O_I12) <- N_age
  dim(O_I32) <- N_age
  dim(O_I42) <- N_age
  dim(O_I13) <- N_age
  dim(O_I23) <- N_age
  dim(O_I43) <- N_age
  dim(O_I14) <- N_age
  dim(O_I24) <- N_age
  dim(O_I34) <- N_age
  dim(O_I231) <- N_age
  dim(O_I241) <- N_age
  dim(O_I341) <- N_age
  dim(O_I132) <- N_age
  dim(O_I142) <- N_age
  dim(O_I342) <- N_age
  dim(O_I123) <- N_age
  dim(O_I143) <- N_age
  dim(O_I243) <- N_age
  dim(O_I124) <- N_age
  dim(O_I134) <- N_age
  dim(O_I234) <- N_age
  dim(O_I2341) <- N_age
  dim(O_I1342) <- N_age
  dim(O_I1243) <- N_age
  dim(O_I1234) <- N_age
  
  dim(O_R1) <- N_age
  dim(O_R2) <- N_age
  dim(O_R3) <- N_age
  dim(O_R4) <- N_age
  dim(O_R12) <- N_age
  dim(O_R13) <- N_age
  dim(O_R14) <- N_age
  dim(O_R23) <- N_age
  dim(O_R24) <- N_age
  dim(O_R34) <- N_age
  dim(O_R123) <- N_age
  dim(O_R124) <- N_age
  dim(O_R134) <- N_age
  dim(O_R234) <- N_age
  dim(O_R1234) <- N_age
  

  
  
  
  
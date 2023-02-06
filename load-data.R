library(tidyverse)
library(lubridate)



df_obs = local({


  source('utils.R')
  # state-level data --------------------------------------------------------

  welfare = openxlsx::read.xlsx(
    xlsxFile = 'https://cpr.uky.edu/sites/ukcpr/files/UKCPR_National_Welfare_Data_Update_022322.xlsx',
    sheet = 'Data'
  ) %>%
    select(state = state_name, year, welf = `AFDC/TANF_FS.4-Person.Benefit`)

  abbreviate_states = function(statename) ifelse(statename == 'District of Columbia', 'DC', state.abb[which(state.name == statename)])

  #these files are available at https://www.americashealthrankings.org/explore/annual as of 7.22.22
  public_health_funding = list.files(path = 'ahr-data') %>%
    map(~read_csv(paste0('ahr-data/', .))) %>%
    map(~mutate(., state = sapply(.$`State Name`, abbreviate_states))) %>%
    map(filter, `Measure Name` == 'Public Health Funding') %>%
    map(select, state, year=Edition, phf=Rank) %>%
    bind_rows()


  minimum_wage = readxl::read_excel('data/publicly-available/mw/mw_state_monthly.xlsx') %>%
    separate(`Monthly Date`, into = c('year','month'), sep='m') %>%
    filter(year>2005) %>%
    select(state = `State Abbreviation`, year, month, fed_mw = `Monthly Federal Minimum`, state_mw = `Monthly State Minimum`) %>%
    mutate(mw = pmax(fed_mw, state_mw))




  policy_data_raw = read_csv('data/publicly-available/covariate_data_from_python.csv')
  regards_cvd_raw = read_csv('data/regards/regards_from_python.csv')[,-1]
  regards_covar_raw = read_csv('data/regards/regards_blcovs_frompython.csv')
  regards_raw = inner_join(regards_cvd_raw, regards_covar_raw)







  regards = regards_raw %>%
    filter(cvdft18dt > as.Date('2007-12-31')) %>%
    make_periodwise(time='cvdft18dt', event='cvdft18', ltfu='ltfu',level='month') %>%
    mutate(white = 1*(Race=='W')) %>%
    group_by(id) %>%
    mutate(lag_event = lag(event)) %>%
    ungroup() %>%
    mutate(hypert = 1*(SBP>=140 | DBP>=90),
           age_08 = Age + (2008 - year(INTDATE)),
           age_08_std = (age_08 - mean(age_08, na.rm=TRUE))/sd(age_08, na.rm=TRUE),
           white = 1*(Race=='W'),
           female = 1*(Gender=='F'),
           cvd = 1*(CAD_SR_ECG == 'Y'),
           smoke = 1*(Smoke_current=='Y'),
           ed0 = 1*(ED_Cat == 'Less than high school'),
           ed1 = 1*(ED_Cat == 'College graduate and above'),
           ed2 = 1*(ED_Cat == 'Some college'),
           srh0 = 1*(Gen_SR_Health == 'Excellent'),
           srh1 = 1*(Gen_SR_Health == 'Very good'),
           srh2 = 1*(Gen_SR_Health == 'Fair'),
           srh3 = 1*(Gen_SR_Health == 'Poor'),
           inc_lt_15 = 1*(Income < 4),
           cvddeath = cvdft18,
           death = death18,
           before_08 = death18dt <= as.Date('2008-01-01')
    ) %>%
    select(id, State, period, event, lag_event, noncvddeath, age_08_std,
           female, white, ed0:ed2, edu=ED_Cat, PSS, hypert, cvd,
           srh0:srh3, srh=Gen_SR_Health, Race, Gender, Age, smoke, Weight,
           REGION, inc_lt_15:before_08) %>%
    set_names(tolower(names(.)))

  regards_bl = regards_raw %>%
    mutate(insamp = 1*(id %in% regards$id))


  policy_data <- policy_data_raw %>%
    mutate(period = 12*(year - 2008) + month - 1) %>% #period=months since Jan 1, 2008
    group_by(state, year, month, period) %>%
    summarise(mw=first(mw), fed_mw=first(Fed_mw),
              welf=first(afdc_snap_4person)/1200,
              phf=first(phunding)) %>%
    ungroup() %>%
    filter(!state %in% c("AK","HI")) %>% #these states aren't in regards
    left_join(filter(., state=='NY') %>% select(period, ny_mw=mw)) %>%
    mutate(mw_g = ifelse(period >= 18 & mw < ny_mw, ny_mw, mw)) %>% #federal=ny mw first in period 18
    group_by(period) %>%
    mutate(mw_scaled = (mw - mean(mw))/sd(mw),
           mw_scaled_g =(mw_g - mean(mw_g))/sd(mw_g)) %>%
    group_by(state) %>%
    mutate(period_std = period / 12,
           a = mw-mw_g, #defining the exposure as the difference between the MW and what it would be under intervention
           a_delagged = a-dplyr::lag(a, n=12), #subtracting 6-month lags to address collinearity w/ time
           a_delagged_per = a_delagged*period_std, #just need to remember that the interactions with period are standardized
           mw_10c = round(10*(mw - fed_mw)),
           mw_bin = 1*(mw_10c==0),
           mw_pois = mw_10c - 1, #always remember to subset the data to mw_bin==0 when modeling this variable
           mw_cummean_lag = lag(cummean(mw_scaled)),
           mw_cummean_lag_g = lag(cummean(mw_scaled_g)),
           a_cummean_lag = mw_cummean_lag - mw_cummean_lag_g) %>%
    ungroup()

  regards %>%
    inner_join(policy_data) %>%
    select(id:smoke, mw, fed_mw, mw_g, period_std:a_cummean_lag, month, year, welf, phf, sweight=weight)

})


df_ipw <- df_obs %>% #annual dataset for the weights (ideally every 6 months since that is almost max frequency of changes, but start here)
  filter(month == 1, period>0) %>% #this is the same dataset used to fit the weights in yearly data
  mutate(mw_qb = qbin(mw, n=10),
         mw_qb_fc = factor(mw_qb, ordered=TRUE),
         mw_10c = round(10*(mw - fed_mw)),
         agecat = qbin(age_08_std, n=5, na.rm=TRUE),
         # demcat = qbin(dem, n=5, na.rm=TRUE),
         # povcat = qbin(pov, n=5, na.rm=TRUE),
         # votcat = qbin(vot, n=5, na.rm=TRUE),
         phfcat = qbin(phf, n=5, na.rm=TRUE),
         welfcat = qbin(welf, n=5, na.rm=TRUE),
         # inccat = qbin(inc, n=5, na.rm=TRUE),
         # ahrcat = qbin(ahr, n=5, na.rm=TRUE),
         # infcat = qbin(infmort, n=5, na.rm=TRUE)
         ) %>%
  select(-pss)

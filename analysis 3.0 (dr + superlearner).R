library(tidyverse)
library(rlang)
library(SuperLearner)
library(rsample)
source('load-data.R')
source('penalized regression helpers.R')
source("../../01 Parallel trends g-formula/IF estimators/superlearner_helpers.R")
source("../../01 Parallel trends g-formula/IF estimators/crossfit estimator.R")
source('sl_aipw_wrappers.R')

intervene = function(df) df %>% mutate(across(starts_with('a'), ~1)) #possibly edit this later to incorporate lagged values, actually MW values etc.

# ipw (used for both aipw and tmle) ----------------------------------------------------------------



df_indiv = df_obs %>%
  left_join(filter(., state=='NY') %>% select(period, ny_mw=mw) %>% unique()) %>%
  arrange(id, period) %>%
  group_by(id) %>%
  mutate(mw_star      =  fed_mw*(period<71) + ny_mw*(period >=71),
         under_regime = cumprod(1*(mw >= mw_star)),#some states go in and out - in this version, we will only use them once.
         under_regime_lag = lag(under_regime, default = 1),
         mw_cummean_lag = lag(cummean(mw)),
         welf1=welf,
         welf2=lag(welf, 12),
         welf3=lag(welf, 24),
         welf4=lag(welf, 36),
         welf5=lag(welf, 48)) %>% #up to five years of lags of welfare
  ungroup()

diverge_periods <- with(df_indiv, unique(period[under_regime_lag== 1 & under_regime == 0])) %>% sort(decreasing = TRUE) %>% c(max(.) + 1, .)

df_ipw_indiv = df_indiv %>%
  filter(period %in% diverge_periods, under_regime_lag==1) %>%
  mutate(period_scaled = (period - mean(period))/sd(period),
         period2_scaled = (period^2  - mean(period^2))/sd(period^2),
         welf_scaled = (welf - mean(welf))/sd(welf)) %>%
  bind_cols(ns(.$age_08_std, df=3) %>%
              as_tibble() %>%
              set_names(paste0('age_ns', names(.))) %>%
              mutate(across(everything(), ~(. - mean(., na.rm=TRUE))/sd(., na.rm=TRUE)))) %>%
  select(id, period, under_regime, under_regime_lag, event, lag_event, female:ed2, starts_with('age_ns'), hypert:srh3, smoke, period_scaled:welf_scaled, sweight)


num_mod = glmnet_f(formula = under_regime ~ white*(period_scaled + period2_scaled),
                   data = df_ipw_indiv,
                   family = 'binomial',
                   weights = sweight,
                   subset = under_regime_lag == 1,
                   alpha = 0,
                   lambda=1/(2 * 1.5^2)) #pre-selected, so no concern about convergence rates

den_mod = glmnet_f(formula = under_regime ~ white*(period_scaled + period2_scaled + welf_scaled + female + ed0+ed1+ed2 + age_ns1+age_ns2+age_ns3 + hypert + cvd + srh1+srh2+srh3 + smoke + welf_scaled),
                   data = df_ipw_indiv,
                   family='binomial',
                   weights = sweight,
                   subset = under_regime_lag == 1,
                   alpha=0,
                   lambda = 1/(2* 1.5^2))


keep_indices = apply(model.frame(den_mod$formula, df_ipw_indiv, na.action = na.pass), 1, function(x) all(!is.na(x)))

df_ipw = df_ipw_indiv %>%
  filter(keep_indices) %>%
  mutate(num_pred = predict(num_mod, newx=model.matrix(num_mod$formula, intervene(.)), type='response')[,1],
         den_pred = predict(den_mod, newx=model.matrix(den_mod$formula, intervene(.)), type='response')[,1]) %>%
  group_by(id) %>%
  mutate(g = cumprod(den_pred),
         ipw = cumprod(num_pred / den_pred)) %>%
  ungroup() %>%
  select(id, period, a=under_regime, under_regime_lag, g, ipw)


# aipw estimator ------------------------------------------------------------

# ice_periods are the periods of latest history of variables in the conditioning event.
# For example, at times 71:82, we are at the outermost expectation already.
# At times 83-107, the outermost expectation conditions on time 82.
ice_periods = diverge_periods[c(-1,-6)] -1
ice_intervals = diff(ice_periods) %>% set_names(ice_periods[-4])
ice_first_lag = c(71:131) %>% set_names() %>% map_dbl(~max(ice_periods[ice_periods<.])) %>% ifelse(is.infinite(.), NA, .)



df_ice = df_indiv %>% filter(period >= 71) %>%
  group_by(period) %>%
  mutate(sweight = sweight/sum(sweight)) %>%
  ungroup() %>%
  mutate(deltay = (event - lag_event + 1)/2) %>%  #transform differences to lie between zero and one
  bind_cols(ns(.$age_08_std, df=3) %>%
              as_tibble() %>%
              set_names(paste0('age_ns', names(.))) %>%
              mutate(across(everything(), ~(. - mean(., na.rm=TRUE))/sd(., na.rm=TRUE)))) %>%
  select(id, t=period, a=under_regime, a_lag=under_regime_lag, deltay, event, lag_event, female:ed2, hypert:srh3, smoke, starts_with('age_ns'), sweight,  welf1:welf5) %>%
  mutate(t_lag1 = ice_first_lag[as.character(t)],
         t_lag2 = t_lag1 + ice_intervals[as.character(t_lag1)],
         t_lag3 = t_lag2 + ice_intervals[as.character(t_lag2)]) %>%
  left_join(select(., id, t, a, welf1:welf5) %>% rename_with(~paste0(.,'_lag1'), t:welf5)) %>%
  left_join(select(., id, t, a, welf1:welf5) %>% rename_with(~paste0(.,'_lag2'), t:welf5)) %>%
  left_join(select(., id, t, a, welf1:welf5) %>% rename_with(~paste0(.,'_lag3'), t:welf5))


qmods_aipw = list()
df_aipw = df_ice %>%
  mutate(Qtt_3 = predict(qmods_aipw[['t3']] <<- glm(deltay ~ white + ns(t-71, df=3) + welf1
                                               + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                               quasibinomial, data=., subset=a==1, weights=sweight),
                         type='response', newdata=.),
         Qtt_2 = predict(qmods_aipw[['t2']] <<- glm(Qtt_3 ~ white + ns(t-t_lag1, df=3) + welf1_lag1 + factor(t_lag1)
                                               + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                               quasibinomial, data=., subset=a_lag1==1, weights=sweight),
                         type='response', newdata=.),
         Qtt_1 = predict(qmods_aipw[['t1']] <<- glm(Qtt_2 ~ white + ns(t-t_lag2, df=3) + welf1_lag2 + factor(t_lag2)
                                               + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                               quasibinomial, data=., subset=a_lag2==1, weights=sweight),
                         type='response', newdata=.),
         Qtt_0 = predict(qmods_aipw[['t0']] <<- glm(Qtt_1 ~ white + ns(t-t_lag3, df=2) + welf1_lag3
                                               + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                               quasibinomial, data=., subset=a_lag3==1, weights=sweight),
                         type='response', newdata=.),
         Qtt = coalesce(Qtt_0, Qtt_1, Qtt_2, Qtt_3)*2 - 1) %>%
  select(id, sweight, t, a, white, deltay, lag_event, starts_with('Q'))



df_eif = df_ipw %>%
  group_by(id) %>%
  mutate(g3=g, g2=lag(g), g1 = lag(g, 2), g0 = lag(g, 3),
         I3=a, I2=lag(a), I1 = lag(a, 2), I0 = lag(a, 3)) %>%
  right_join(df_aipw %>% rename(period=t)) %>%
  arrange(id, period) %>%
  select(id, period, white, a, deltay,
         ylag=lag_event, starts_with('Q'), starts_with('g'), starts_with('I')) %>%
  fill(g:I0) %>%
  ungroup() %>%
  mutate(
    #replacing NAs for the Q and I() functions with 0 and for the g functions with 1
    #allows these to only enter the phi_tt and phi_ttmin1 expressions below
    #when they are actually needed. So for example, for period <= 83, we only
    #need (y-Qtt_3)/g + Qtt.
    across(starts_with('Q'), replace_na, 0),
    across(starts_with('g'), replace_na, 1),
    across(starts_with('I'), replace_na, 0),
    phi_tt = (deltay - Qtt_3)*I3/g3 + (Qtt_3 - Qtt_2)*I2/g2 + (Qtt_2 - Qtt_1)*I1/g1 + (Qtt_1 - Qtt_0)*I0/g0 + Qtt) %>%
  left_join(filter(., period==71) %>% select(id, phi_00 = ylag))

#all the below needs to take into account sampling weights
eif_obs_rates = df_obs %>%
  group_by(white, period) %>%
  mutate(rhat = sum(event)/n(),
         eif_rhat = event - rhat) %>%
  ungroup() %>%
  select(id, white,period,rhat,eif_rhat)

eif_cf_rates = df_eif %>%
  group_by(id) %>%
  mutate(psihat = phi_00 + cumsum(phi_tt)) %>%
  group_by(period, white) %>%
  mutate(rhat     = mean(psihat),
         eif_rhat = psihat - rhat) %>%
  ungroup() %>%
  select(id, period, white, rhat, eif_rhat)



df_risk_aipw = list(observed = eif_obs_rates %>%  filter(period < 71),
                   counterfactual = eif_cf_rates) %>%
  bind_rows(.id='world') %>%
  arrange(id, white, period) %>%
  group_by(id) %>%
  mutate(Rhat = cumsum(rhat), #rhat is rate, Rhat is risk
         eif_Rhat = cumsum(eif_rhat)) %>%
  group_by(white, period) %>%
  summarise(est_Rhat = mean(Rhat),
            var_Rhat = mean(eif_Rhat^2)/n(),
            lwr = est_Rhat - 1.96*sqrt(var_Rhat), upr = est_Rhat + 1.96*sqrt(var_Rhat))

ggplot(df_risk_aipw, aes(x=period, y=est_Rhat, color=factor(white), fill=factor(white))) +
  geom_step() +
  pammtools::geom_stepribbon(aes(ymin=lwr, ymax=upr), alpha=0.2, color=NA) +
  coord_cartesian(ylim=c(0,0.12))




# aipw w/ superlearner ----------------------------------------------------



qmods_aipw <- list()
q_library = c('SL.mean', 'SL.glmnet')

a=predict(qmods_aipw[['Q2']] <- sl_wrapper(deltay ~ .,  data = df_ice %>% select(t:sweight, t_lag1),
             subset=a_lag==1,  SL.library = q_library,  weights = sweight,
             cvControl=list(V=4)), newdata=df_ice, type='response')

qmods_aipw[['Q3']] <-sl_wrapper(deltay ~ .,  data = df_ice %>% select(t:sweight),
                                 subset=a==1,  SL.library = q_library,  weights = sweight,
                                 cvControl=list(V=4))



#step 1. one Qt at a time (sampled data set)
df_test = df_ice[sample(1:nrow(df_ice), size=1e4), ] %>% select(id, event, lag_event, female, white, sweight, starts_with('a'), starts_with('welf'), starts_with('t'))# %>% mutate(event = rbinom(n(), size=1, prob = plogis(-2 + a + female + white)))
qmods_aipw[['Q3']] <-sl_wrapper_aipw(event ~ ., data=df_test, family=binomial, event, welf1:welf5, a)
df_test$Q3t = predict_aipw(qmods_aipw[['Q3']], data=df_test)
qmods_aipw[['Q2']] <-sl_wrapper_aipw(Q3t ~ ., data=df_test, family=gaussian(link='logit'), Q3t, welf1_lag1:welf5_lag1, a_lag1, t_lag1)
df_test$Q2t = predict_aipw(qmods_aipw[['Q2']], data=df_test)
qmods_aipw[['Q1']] <-sl_wrapper_aipw(Q2t ~ ., data=df_test, family=gaussian(link='logit'), Q2t, welf1_lag2:welf5_lag2, a_lag2, t_lag2)
df_test$Q1t = predict_aipw(qmods_aipw[['Q1']], data=df_test)
qmods_aipw[['Q0']] <-sl_wrapper_aipw(Q1t ~ ., data=df_test, family=gaussian(link='logit'), Q1t, welf1_lag3:welf5_lag3, a_lag3, t_lag3)
df_test$Q0t = predict_aipw(qmods_aipw[['Q0']], data=df_test)

qmods_aipw[['Q3tmin1']] <-sl_wrapper_aipw(lag_event ~ ., data=df_test, family=binomial, lag_event, welf1:welf5, a)
df_test$Q3tmin1 = predict_aipw(qmods_aipw[['Q3']], data=df_test)
qmods_aipw[['Q2tmin1']] <-sl_wrapper_aipw(Q3tmin1 ~ ., data=df_test, family=gaussian(link='logit'), Q3tmin1, welf1_lag1:welf5_lag1, a_lag1, t_lag1)
df_test$Q2tmin1 = predict_aipw(qmods_aipw[['Q2']], data=df_test)
qmods_aipw[['Q1tmin1']] <-sl_wrapper_aipw(Q2tmin1 ~ ., data=df_test, family=gaussian(link='logit'), Q2tmin1, welf1_lag2:welf5_lag2, a_lag2, t_lag2)
df_test$Q1tmin1 = predict_aipw(qmods_aipw[['Q1']], data=df_test)
qmods_aipw[['Q0tmin1']] <-sl_wrapper_aipw(Q1tmin1 ~ ., data=df_test, family=gaussian(link='logit'), Q1tmin1, welf1_lag3:welf5_lag3, a_lag3, t_lag3)
df_test$Q0tmin1 = predict_aipw(qmods_aipw[['Q0']], data=df_test)



#step 2: all at once using pipes (sampled data set)
df_aipw_sl = df_test %>%
  mutate(Q3t = predict_aipw(qmods_aipw[['Q3t']] <<- sl_wrapper_aipw(event ~ ., data=., family=binomial, library=q_library, event, welf1:welf5, a),
                          data = .)) %>%
  mutate(Q3tmin1 = predict_aipw(qmods_aipw[['Q3tmin1']] <<- sl_wrapper_aipw(lag_event ~ ., data=., family=binomial, library=q_library, event, welf1:welf5, a),
                                data = .),
         Q3 = (Q3t - Q3tmin1 + 1)/2) %>%
  mutate(Q2 = predict_aipw(qmods_aipw[['Q2']] <<- sl_wrapper_aipw(Q3 ~., data=., family=gaussian(link='logit'), library=q_library, Q3, welf1_lag1:welf5_lag1, a_lag1, t_lag1),
                            data=.)) %>%
  mutate(Q1 = predict_aipw(qmods_aipw[['Q1']] <<- sl_wrapper_aipw(Q2 ~., data=., family=gaussian(link='logit'), library=q_library, Q2, welf1_lag2:welf5_lag2, a_lag2, t_lag2),
                            data=.)) %>%
  mutate(Q0 = predict_aipw(qmods_aipw[['Q0']] <<- sl_wrapper_aipw(Q1~., data=., family=gaussian(link='logit'), library=q_library, Q1, welf1_lag3:welf5_lag3, a_lag3, t_lag3),
                            data=.),
         Q = coalesce(Q0,Q1,Q2,Q3)*2 - 1) %>%
  select(id, sweight, t, a, white, event, lag_event, starts_with('Q'))


#step 3: full data set
df_aipw_sl = df_ice %>%
  mutate(Q3t = predict_aipw(qmods_aipw[['Q3t']] <<- sl_wrapper_aipw(event ~ ., data=., family=binomial, library=q_library, event, welf1:welf5, a),
                            data = .)) %>%
  mutate(Q2t = predict_aipw(qmods_aipw[['Q2t']] <<- sl_wrapper_aipw(Q3t ~., data=., family=gaussian(link='logit'), library=q_library, Q3t, welf1_lag1:welf5_lag1, a_lag1, t_lag1),
                            data=.)) %>%
  mutate(Q1t = predict_aipw(qmods_aipw[['Q1t']] <<- sl_wrapper_aipw(Q2t ~., data=., family=gaussian(link='logit'), library=q_library, Q2t, welf1_lag2:welf5_lag2, a_lag2, t_lag2),
                            data=.)) %>%
  mutate(Q0t = predict_aipw(qmods_aipw[['Q0t']] <<- sl_wrapper_aipw(Q1t~., data=., family=gaussian(link='logit'), library=q_library, Q1t, welf1_lag3:welf5_lag3, a_lag3, t_lag3),
                            data=.),
         Qt = coalesce(Q0t,Q1t,Q2t,Q3t)) %>%
  mutate(Q3tmin1 = predict_aipw(qmods_aipw[['Q3tmin1']] <<- sl_wrapper_aipw(lag_event ~ ., data=., family=binomial, library=q_library, event, welf1:welf5, a),
                                data = .)) %>%
  mutate(Q2tmin1 = predict_aipw(qmods_aipw[['Q2tmin1']] <<- sl_wrapper_aipw(Q3tmin1 ~., data=., family=gaussian(link='logit'), library=q_library, Q3t, welf1_lag1:welf5_lag1, a_lag1, t_lag1),
                                data=.)) %>%
  mutate(Q1tmin1 = predict_aipw(qmods_aipw[['Q1tmin1']] <<- sl_wrapper_aipw(Q2tmin1 ~., data=., family=gaussian(link='logit'), library=q_library, Q2t, welf1_lag2:welf5_lag2, a_lag2, t_lag2),
                                data=.)) %>%
  mutate(Q0tmin1 = predict_aipw(qmods_aipw[['Q0tmin1']] <<- sl_wrapper_aipw(Q1tmin1~., data=., family=gaussian(link='logit'), library=q_library, Q1t, welf1_lag3:welf5_lag3, a_lag3, t_lag3),
                                data=.),
         Qtmin1 = coalesce(Q0tmin1,Q1tmin1,Q2tmin1,Q3tmin1)) %>%
  select(id, sweight, t, a, white, event, lag_event, starts_with('Q'))


#step 4: if this works delete step 3
df_aipw_sl = df_ice %>%
  mutate(Q3t = predict_aipw(qmods_aipw[['Q3t']] <<- sl_wrapper_aipw(event ~ ., data=., family=binomial, event, welf1:welf5, a),
                            data = .),
         Q3tmin1 = predict_aipw(qmods_aipw[['Q3tmin1']] <<- sl_wrapper_aipw(lag_event ~ ., data=., family=binomial, lag_event, welf1:welf5, a),
                                data = .),
         Q3 = (Q3t - Q3tmin1 + 1)/2) %>%
  mutate(Q2 = predict_aipw(qmods_aipw[['Q2']] <<- sl_wrapper_aipw(Q3 ~., data=., family=gaussian(link='logit'), Q3, welf1_lag1:welf5_lag1, a_lag1, t_lag1),
                         data=.)) %>%
  mutate(Q1 = predict_aipw(qmods_aipw[['Q1']] <<- sl_wrapper_aipw(Q2 ~., data=., family=gaussian(link='logit'), Q2, welf1_lag2:welf5_lag2, a_lag2, t_lag2),
                           data=.)) %>%
  mutate(Q0 = predict_aipw(qmods_aipw[['Q0']] <<- sl_wrapper_aipw(Q1~., data=., family=gaussian(link='logit'), Q1, welf1_lag3:welf5_lag3, a_lag3, t_lag3),
                           data=.),
         Q = coalesce(Q0,Q1,Q2,Q3)) %>%
  select(id, sweight, t, a, white, event, lag_event, starts_with('Q'))

#step 5: package for cross fitting
nuissance = function(q_library, g_library) {
  list(
    Q3t = function(df) sl_wrapper_aipw(event ~ .,    data=df, family=binomial(), library=q_library, event, welf1:welf5, a),
    Q3tmin1 = function(df) sl_wrapper_aipw(lag_event ~ ., data=df, family=binomial(), library=q_library, lag_event, welf1:welf5, a),
    Q3 = function(df) pseudomodel( (Q3t - Q3tmin1 + 1)/2, data=df), #note scale change
    Q2 = function(df) sl_wrapper_aipw(Q3 ~ ., data=df, family=gaussian('logit'), library=q_library, Q3, welf1_lag1:welf5_lag1, a_lag1, t_lag1),
    Q1 = function(df) sl_wrapper_aipw(Q2 ~ ., data=df, family=gaussian('logit'), library=q_library, Q2, welf1_lag2:welf5_lag2, a_lag2, t_lag2),
    Q0 = function(df) sl_wrapper_aipw(Q1~ .,  data=df, family=gaussian('logit'), library=q_library, Q1, welf1_lag3:welf5_lag3, a_lag3, t_lag3),
    g = function(df) sl_wrapper_aipw(a ~ ., data=df, family=binomial(), library=g_library, welf1:welf5, a, a_lag1)
  )
}

fit_nuiss = function(df_cv, nuissance) {
  df_cv %>%
    mutate(nuissance_fits = map(splits, ~get_nuissance_fits(data=analysis(.),
                                                            nuissance_models = nuissance,
                                                            first_colname = 'y2')))
}

nuissance_models = nuissance(q_library=c('SL.glm'), g_library=c('SL.glm'))

#step 5a: one model at a time
qmods_aipw[['Q3t']] <- nuissance_models$Q3t(df_test)
qmods_aipw[['Q3tmin1']] <- nuissance_models$Q3tmin1(df_test)
df_test$Q3t = predict_aipw(qmods_aipw$Q3t, df_test)
df_test$Q3tmin1 = predict_aipw(qmods_aipw$Q3tmin1, df_test)
qmods_aipw[['Q3']] <- nuissance_models$Q3(df_test)
df_test$Q3 = predict_aipw(qmods_aipw$Q3, data=df_test)
qmods_aipw[['Q2']] <-nuissance_models$Q2(df_test)
df_test$Q2 = predict_aipw(qmods_aipw[['Q2']], data=df_test)
qmods_aipw[['Q1']] <- nuissance_models$Q1(df_test)
df_test$Q1 = predict_aipw(qmods_aipw[['Q1']], data=df_test)
qmods_aipw[['Q0']] <-nuissance_models$Q0(df_test)
df_test$Q0t = predict_aipw(qmods_aipw[['Q0']], data=df_test)

#step 5b: all models using meta-function
nuissance_fits = get_nuissance_fits(data = df_test, nuissance_models = nuissance_models, first_colname = 'event', return_data = TRUE)

#step 6: cross fitting
n_folds = 2


df_aipw_cv = df_ice %>%
  group_vfold_cv(group=id, v=n_folds) %>%
  mutate(nuissance_fits = map(splits,
                              ~get_nuissance_fits(data=analysis(.),
                                                  nuissance_models = nuissance_models,
                                                  first_colname = 'event')))

calc_eif_pooled = function(df) {
  df %>%
    group_by(id) %>%
    arrange(id, t) %>%
    mutate(g = cumprod(g)) %>%
    left_join(select(., id, t_lag1 = t, g_lag1 = g)) %>%
    left_join(select(., id, t_lag2 = t, g_lag2 = g)) %>%
    left_join(select(., id, t_lag3 = t, g_lag3 = g)) %>%
    left_join(select(filter(., t==71), id, event_0 = lag_event)) %>%
    mutate(Q2 = coalesce(Q2, Q3),
           Q1 = coalesce(Q1, Q2, Q3),
           Q0 = coalesce(Q0, Q1, Q2, Q3),
           across(starts_with('Q'), ~(.*2 - 1)),
           across(starts_with('g_lag'), ~ifelse(is.na(.), 1, .)),
           across(starts_with('a_lag'), ~ifelse(is.na(.), 0, .)), #this causes all irrelevant Q's to zero out of the influence function
           phi_t = (1*(a==1)/g)*(event - lag_event - Q3)
           + (1*(a_lag1==1)/g_lag1)*(Q3 - Q2)
           + (1*(a_lag2==1)/g_lag2)*(Q2 - Q1)
           + (1*(a_lag3==1)/g_lag3)*(Q1 - Q0) + Q0,
           psi_t = event_0 + cumsum(phi_t)
    )
}

df_aipw_cv = df_aipw_cv %>%
  mutate(outsamp_nuissance_estimates = map2(nuissance_fits, splits,
                                            ~reduce2(.x = .x,
                                                     .y = parse_exprs(names(.x)),
                                                     .f = add_outsamp_pred,
                                                     intervene = intervene,
                                                     .init = assessment(.y))))

df_eif = df_aipw_cv %>%
  select(fold=id, outsamp_nuissance_estimates) %>%
  unnest(outsamp_nuissance_estimates) %>%
  calc_eif_pooled()


df_{left_join(summarise_psi(., .fns=list(est=mean, var=var)), zivich_estimator(.))} # do the variance both ways.



qmods_aipw_cf = as.list(rep(NA, nrow(df_aipw_cf))) %>% set_names(df_aipw_cf$id) #one idea



df_aipw %>%
  mutate(Qtt_3 = predict(qmods_aipw[['t3']] <<- glm(deltay ~ white + ns(period-71, df=3) + welf_scaled
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a==1, weights=sweight),
                         type='response', newdata=.),
         Qtt_2 = predict(qmods_aipw[['t2']] <<- glm(Qtt_3 ~ white + ns(period-t_lag1, df=3) + welf_lag1 + factor(t_lag1)
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag1==1, weights=sweight),
                         type='response', newdata=.),
         Qtt_1 = predict(qmods_aipw[['t1']] <<- glm(Qtt_2 ~ white + ns(period-t_lag2, df=3) + welf_lag2 + factor(t_lag2)
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag2==1, weights=sweight),
                         type='response', newdata=.),
         Qtt_0 = predict(qmods_aipw[['t0']] <<- glm(Qtt_1 ~ white + ns(period-t_lag3, df=2) + welf_lag3
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag3==1, weights=sweight),
                         type='response', newdata=.),
         Qtt = coalesce(Qtt_0, Qtt_1, Qtt_2, Qtt_3)*2 - 1) %>%
  select(id, sweight, period, a, white, deltay, lag_event, starts_with('Q'))



# tmle 1 - separate models for t and t-1 --------------------------------------------------------------------

# TMLE is still predicting rates outside of the range. One idea is to take
# event - lag_event and transform it to lie between 0 and 1 and then estimate
# a quasibinomial model, instead of all the separate models. So,
# (event - lag_event) lies between -1 and 1, so add 1 and divide by 2.


qmods_tmle = list()
df_tmle = df_ipw %>%
  group_by(id) %>%
  mutate(g3=g, g2=lag(g), g1 = lag(g, 2), g0 = lag(g, 3),
         I3=a, I2=lag(a), I1 = lag(a, 2), I0 = lag(a, 3)) %>%
  right_join(df_ice) %>%
  arrange(id, period) %>%
  fill(g:I0) %>%
  ungroup() %>%
  mutate(
    deltay = (event - lag_event + 1)/2, #transform differences to lie between zero and one
    #replacing NAs for the I() functions with 0 and for the g functions with 1, which
    # just makes it transparent that individuals with NA for those don't enter the model
    across(starts_with('g'), replace_na, 1),
    across(starts_with('I'), replace_na, 0)) %>%
  mutate(Qtt_3 = predict(qmods_tmle[['t3']] <<- glm(event ~ white + ns(period-71, df=3) + welf_scaled
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_3 = predict(qmods_tmle[['t3_eps']] <<- glm(event ~ white + offset(plogis(Qtt_3)), data=., weights = sweight * I3/g3 ),
                         type='response', newdata=.),
         Qtt_2 = predict(qmods_tmle[['t2']] <<- glm(Qtt_3 ~ white + ns(period-t_lag1, df=3) + welf_lag1 + factor(t_lag1)
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag1==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_2 = predict(qmods_tmle[['t2_eps']] <<- glm(Qtt_3 ~ white + offset(plogis(Qtt_2)), data=., weights = sweight * I2/g2 ),
                         type='response', newdata=.),
         Qtt_1 = predict(qmods_tmle[['t1']] <<- glm(Qtt_2 ~ white + ns(period-t_lag2, df=3) + welf_lag2 + factor(t_lag2)
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag2==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_1 = predict(qmods_tmle[['t1_eps']] <<- glm(Qtt_2 ~ white + offset(plogis(Qtt_1)), data=., weights = sweight * I1/g1),
                         type='response', newdata=.),
         Qtt_0 = predict(qmods_tmle[['t0']] <<- glm(Qtt_1 ~ white + ns(period-t_lag3, df=2) + welf_lag3
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag3==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_0 = predict(qmods_tmle[['t0_eps']] <<- glm(Qtt_1 ~ white + offset(plogis(Qtt_0)), data=., weights = sweight * I0/g0),
                         type='response', newdata=.),
         Qtt = coalesce(Qtt_0, Qtt_1, Qtt_2, Qtt_3),
         Qttmin1_3 = predict(qmods_tmle[['t3min1']] <<- glm(event ~ white + ns(period-71, df=3) + welf_scaled
                                                            + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                            quasibinomial, data=., subset=a==1, weights=sweight),
                             type='response', newdata=.)) %>%
  mutate(Qttmin1_3 = predict(qmods_tmle[['t3_eps']] <<- glm(lag_event ~ white + offset(plogis(Qttmin1_3)), data=., weights = sweight * I3/g3 ),
                             type='response', newdata=.),
         Qttmin1_2 = predict(qmods_tmle[['t2min1']] <<- glm(Qttmin1_3 ~ white + ns(period-t_lag1, df=3) + welf_lag1 + factor(t_lag1)
                                                            + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                            quasibinomial, data=., subset=a_lag1==1, weights=sweight),
                             type='response', newdata=.)) %>%
  mutate(Qttmin1_2 = predict(qmods_tmle[['t2min1_eps']] <<- glm(Qttmin1_3 ~ white + offset(plogis(Qttmin1_2)), data=., weights = sweight * I2/g2 ),
                             type='response', newdata=.),
         Qttmin1_1 = predict(qmods_tmle[['t1min1']] <<- glm(Qttmin1_2 ~ white + ns(period-t_lag2, df=3) + welf_lag2 + factor(t_lag2)
                                                            + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                            quasibinomial, data=., subset=a_lag2==1, weights=sweight),
                             type='response', newdata=.)) %>%
  mutate(Qttmin1_1 = predict(qmods_tmle[['t1min1_eps']] <<- glm(Qttmin1_2 ~ white + offset(plogis(Qttmin1_1)), data=., weights = sweight * I1/g1),
                             type='response', newdata=.),
         Qttmin1_0 = predict(qmods_tmle[['t0min1']] <<- glm(Qttmin1_1 ~ white + ns(period-t_lag3, df=2) + welf_lag3
                                                            + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                            quasibinomial, data=., subset=a_lag3==1, weights=sweight),
                             type='response', newdata=.)) %>%
  mutate(Qttmin1_0 = predict(qmods_tmle[['t0min1_eps']] <<- glm(Qttmin1_1 ~ white + offset(plogis(Qttmin1_0)), data=., weights = sweight * I0/g0),
                             type='response', newdata=.),
         Qttmin1 = coalesce(Qttmin1_0, Qttmin1_1, Qttmin1_2, Qttmin1_3))



phi00_tmle = df_tmle %>%
  filter(period==71) %>%
  group_by(white) %>%
  summarise(phi00 = mean(lag_event))

#need to take into account sampling weights here
df_psi_tmle = df_tmle %>%
  group_by(white, period) %>%
  summarise(phi_tt = mean(Qtt, na.rm=TRUE),
            phi_ttmin1 = mean(Qttmin1, na.rm=TRUE)) %>%
  left_join(phi00_tmle) %>%
  group_by(white) %>%
  mutate(rhat = phi00 + cumsum(phi_tt - phi_ttmin1))



df_risk_tmle = list(observed = eif_obs_rates %>%
                      filter(period < 71) %>%
                      group_by(white, period) %>%
                      summarise(rhat = mean(rhat)),
                   counterfactual = df_psi_tmle %>%
                     select(white, period, rhat)) %>%
  bind_rows(.id='world') %>%
  arrange(white, period) %>%
  mutate(est_Rhat = cumsum(rhat)) %>%
  left_join(df_risk_aipw %>% select(white, period, var_Rhat)) %>%
  mutate(lwr = est_Rhat - 1.96*sqrt(var_Rhat), upr = est_Rhat + 1.96*sqrt(var_Rhat))


ggplot(df_risk_tmle, aes(x=period, y=est_Rhat, color=factor(white), fill=factor(white))) +
  geom_step() +
  pammtools::geom_stepribbon(aes(ymin=lwr, ymax=upr), alpha=0.2, color=NA) +
  coord_cartesian(ylim=c(0,0.12))


# tmle 2 - model differences ----------------------------------------------

qmods_tmle = list()
df_tmle = df_ipw %>%
  group_by(id) %>%
  mutate(g3=g, g2=lag(g), g1 = lag(g, 2), g0 = lag(g, 3),
         I3=a, I2=lag(a), I1 = lag(a, 2), I0 = lag(a, 3)) %>%
  right_join(df_ice) %>%
  arrange(id, period) %>%
  fill(g:I0) %>%
  ungroup() %>%
  mutate(
    deltay = (event - lag_event + 1)/2, #transform differences to lie between zero and one
    #replacing NAs for the I() functions with 0 and for the g functions with 1, which
    # just makes it transparent that individuals with NA for those don't enter the model
    across(starts_with('g'), replace_na, 1),
    across(starts_with('I'), replace_na, 0)) %>%
  mutate(Qtt_3 = predict(qmods_tmle[['t3']] <<- glm(deltay ~ white + ns(period-71, df=3) + welf_scaled
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_3 = predict(qmods_tmle[['t3_eps']] <<- glm(deltay ~ white + offset(plogis(Qtt_3)), data=., weights = sweight * I3/g3 ),
                         type='response', newdata=.),
         Qtt_2 = predict(qmods_tmle[['t2']] <<- glm(Qtt_3 ~ white + ns(period-t_lag1, df=3) + welf_lag1 + factor(t_lag1)
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag1==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_2 = predict(qmods_tmle[['t2_eps']] <<- glm(Qtt_3 ~ white + offset(plogis(Qtt_2)), data=., weights = sweight * I2/g2 ),
                         type='response', newdata=.),
         Qtt_1 = predict(qmods_tmle[['t1']] <<- glm(Qtt_2 ~ white + ns(period-t_lag2, df=3) + welf_lag2 + factor(t_lag2)
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag2==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_1 = predict(qmods_tmle[['t1_eps']] <<- glm(Qtt_2 ~ white + offset(plogis(Qtt_1)), data=., weights = sweight * I1/g1),
                         type='response', newdata=.),
         Qtt_0 = predict(qmods_tmle[['t0']] <<- glm(Qtt_1 ~ white + ns(period-t_lag3, df=2) + welf_lag3
                                                    + age_ns1+age_ns2+age_ns3+female+hypert+cvd+smoke+srh0+srh1+srh2+srh3+ed0+ed1+ed2,
                                                    quasibinomial, data=., subset=a_lag3==1, weights=sweight),
                         type='response', newdata=.)) %>%
  mutate(Qtt_0 = predict(qmods_tmle[['t0_eps']] <<- glm(Qtt_1 ~ white + offset(plogis(Qtt_0)), data=., weights = sweight * I0/g0),
                         type='response', newdata=.),
         Qtt = coalesce(Qtt_0, Qtt_1, Qtt_2, Qtt_3)*2 - 1)



phi00_tmle = df_tmle %>%
  filter(period==71) %>%
  group_by(white) %>%
  summarise(phi00 = mean(lag_event))

#need to take into account sampling weights here
df_psi_tmle = df_tmle %>%
  group_by(white, period) %>%
  summarise(phi_tt = mean(Qtt, na.rm=TRUE)) %>%
  left_join(phi00_tmle) %>%
  group_by(white) %>%
  mutate(rhat = phi00 + cumsum(phi_tt))



df_risk_tmle = list(observed = eif_obs_rates %>%
                      filter(period < 71) %>%
                      group_by(white, period) %>%
                      summarise(rhat = mean(rhat)),
                    counterfactual = df_psi_tmle %>%
                      select(white, period, rhat)) %>%
  bind_rows(.id='world') %>%
  arrange(white, period) %>%
  mutate(est_Rhat = cumsum(rhat)) %>%
  left_join(df_risk_aipw %>% select(white, period, var_Rhat)) %>%
  mutate(lwr = est_Rhat - 1.96*sqrt(var_Rhat), upr = est_Rhat + 1.96*sqrt(var_Rhat))


ggplot(df_risk_tmle, aes(x=period, y=est_Rhat, color=factor(white), fill=factor(white))) +
  geom_step() +
  pammtools::geom_stepribbon(aes(ymin=lwr, ymax=upr), alpha=0.2, color=NA) +
  coord_cartesian(ylim=c(0,0.12))

# Simulate data to get power analysis #

library(here)
library(data.table)
library(broom)

# set seed 
set.seed(20241123)

load(here("data", "rsv_inselspital_sim.RData"))

# number of repetitions we will do for analyses (aim for 1000 or more)
n.rep <- 1000
plot.rep <- 100
sign.val <- 0.01
plot.year <- 2024

# based on https://nickch-k.github.io/EconometricsSlides/Week_08/Power_Simulations.html

# use historical data to get for each week the percentage that were hospitalised

# data preparation ####
## av. annual birth rate in Kt Bern 2008 - 2019 for the three age groups
pop.kt.bern <- data.table(
  age.grp = c("infant", "toddler", "child"),
  av.birthpop = c(rep(9691, 2), 9691 * 3)
)

## calculate the average hosp.rate per season ####
## using averages over data from 2008 - 2019 (We have admitted notably more patients since 2008) 
### for major season
hosp.rate.maj <- rsv.hist[, .N,  by = .(age.grp, season, season.nr)] |> 
  # calculate average numbers in major seasons
  _[season == "maj", .(N = round(sum(N) / length(unique(season.nr)))),
    by = age.grp] |> 
  # join data with the average population
  _[pop.kt.bern, on = c("age.grp")] |> 
  # calculate average hospitalisation rate per age group and season
  _[, rate := N / av.birthpop]  

### for minor seasons
hosp.rate.min <- rsv.hist[, .N,  by = .(age.grp, season, season.nr)] |> 
  _[season == "min", .(N = round(sum(N) / length(unique(season.nr)))),
    by = age.grp]  |> 
  _[pop.kt.bern, on = c("age.grp")] |> 
  _[, rate := N / av.birthpop]  # av. annual birth rate in Kt Bern 2008 - 2019

## calculate the average hospitalisation rate per week ####
prop.hosp <- rsv.hist[, .N, by = .(adm.isow, age.grp, season, season.nr)]
## using averages over data from 2008 - 2019 (We have admitted notably more patients since 2008)
## calculate the number of cases by week separately in major and minor seasons
prop.hosp.maj <- prop.hosp[season == "maj", ]
prop.hosp.min <- prop.hosp[season == "min", ]
## calculate the total number of admissions per season
prop.hosp.maj[, tot.hosp := sum(N), by = .(age.grp, season.nr)]
prop.hosp.min[, tot.hosp := sum(N), by = .(age.grp, season.nr)]
## for each week calculate how many patients of the total admission number for that season were admitted
prop.hosp.maj[, perc.hosp := N / tot.hosp]
prop.hosp.min[, perc.hosp := N / tot.hosp]


# function to prepare datasets we'll use for comparison ####
# generate data
gen_rsv_data <- function(
    birth.pop = 10000,
    vac.cov = 0,
    vac.eff = 0.8,
    scenario = "maj",
    nirsevimab.avail = "no"
  ) {
  # select which scenario we will use
  if(scenario == "maj") {
    prop.hosp <- prop.hosp.maj
    hosp.rate <- hosp.rate.maj
  } else if(scenario == "min") {
    prop.hosp <- prop.hosp.min
    hosp.rate <- hosp.rate.min
  } else {
    stop("wrong scenario specification, it should be maj or min")
  }
  
  # hospitalisation likelihood if you had nirsevimab
  hosp.rate.infant.nirse <- hosp.rate[age.grp == "infant", rate] * (1 - vac.eff)
  # hospitalisation likelihood if you did not have nirsevimab for all age groups
  hosp.rate.infant <- hosp.rate[age.grp == "infant", rate]
  hosp.rate.toddler <- hosp.rate[age.grp == "toddler", rate]
  hosp.rate.child <- hosp.rate[age.grp == "child", rate]
  # randomly select which distribution we will take
  rand.adm.set <- sample(unique(prop.hosp$season.nr), 1L)
  # select distribution of weekly hospitalisation rate per season
  week.of.admission <- prop.hosp[season.nr == rand.adm.set, ]
  # generate data
  rsv <- data.table(
    # add variable with 'all' children less 5 years
    age.grp = c(
      rep("infant", birth.pop), rep("toddler", birth.pop),
      rep("child", birth.pop * 3L)
    ),
    # add a variable indicating if a child received Nirsevimab or not
    nirsevimab = c(
      # only infants qualify for Nirsevimab
      sample(
        c("yes", "no"), birth.pop, replace = TRUE, prob = c(vac.cov, 1-vac.cov)
      ),
      # older children don't qualify
      rep("no", birth.pop * 4)
    )
  )
  # Add variable indicating if child was hospitalised or not hospitalisation
  # infant without nirsevimab
  rsv[
    age.grp == "infant" & nirsevimab == "no",
    hosp := sample(
      c("yes", "no"), rsv[age.grp == "infant" & nirsevimab == "no", .N],
      replace = TRUE, prob = c(hosp.rate.infant, 1 - hosp.rate.infant)
    )
  ]
  # infant with nirsevimab
  rsv[
    age.grp == "infant" & nirsevimab == "yes",
    hosp := sample(
      c("yes", "no"), rsv[age.grp == "infant" & nirsevimab == "yes", .N],
      replace = TRUE, prob = c(hosp.rate.infant.nirse, 1 - hosp.rate.infant.nirse)  # vaccine efficacy will reduce likelhood of being hospitalised
    )
  ]
  # toddler, does not qualify for nirsevimab
  rsv[
    age.grp == "toddler",
    hosp := sample(
      c("yes", "no"), rsv[age.grp == "toddler", .N],
      replace = TRUE, prob = c(hosp.rate.toddler, 1 - hosp.rate.toddler)
    )
  ]
  # children, do not qualify for nirsevimab
  rsv[
    age.grp == "child",
    hosp := sample(
      c("yes", "no"), rsv[age.grp == "child", .N],
      replace = TRUE, prob = c(hosp.rate.child, 1 - hosp.rate.child)
    )
  ]
  # Add week when were children hospitalised (based on historical distribution)
  rsv[
    age.grp == "infant" & hosp == "yes" & nirsevimab == "no",
    adm.week := sample(
      week.of.admission[age.grp == "infant", adm.isow],
      rsv[age.grp == "infant" & hosp == "yes" & nirsevimab == "no", .N],
      replace = TRUE,
      prob = week.of.admission[age.grp == "infant", perc.hosp]
    )
  ]
  rsv[
    age.grp == "infant" & hosp == "yes" & nirsevimab == "yes",
    adm.week := sample(
      week.of.admission[age.grp == "infant", adm.isow],
      rsv[age.grp == "infant" & hosp == "yes" & nirsevimab == "yes", .N],
      replace = TRUE,
      prob = week.of.admission[age.grp == "infant", perc.hosp]
    )
  ]
  # toddler
  rsv[
    age.grp == "toddler" & hosp == "yes",
    adm.week := sample(
      week.of.admission[age.grp == "toddler", adm.isow],
      rsv[age.grp == "toddler" & hosp == "yes", .N],
      replace = TRUE,
      prob = week.of.admission[age.grp == "toddler", perc.hosp]
    )
  ]
  rsv[
    age.grp == "child" & hosp == "yes",
    adm.week := sample(
      week.of.admission[age.grp == "child", adm.isow],
      rsv[age.grp == "child" & hosp == "yes", .N],
      replace = TRUE,
      prob = week.of.admission[age.grp == "child", perc.hosp]
    )
  ]
  # only return hospitalised patients (simulating hospital data)
  data <- rsv[hosp == "yes", .N, by = .(age.grp, adm.week, nirsevimab)]
  setnames(data, "N", "hosp")
  data[, nirsevimab.avail := nirsevimab.avail]
  
  return(data)
}

prep_analysis_data <- function(birth.pop = 10000,
                               vac.cov = 0,
                               vac.eff = 0.8,
                               scenario = "maj") {
  # first generate 'historic' data when we did not have Nirsevimab
  gen.rsv.hist <- gen_rsv_data(
    birth.pop,
    vac.cov = 0,
    vac.eff,
    scenario,
    nirsevimab.avail = "no"
  )
  # generate 'current' data when Nirsevimab is available
  gen.rsv <- gen_rsv_data(
    birth.pop,
    vac.cov,
    vac.eff,
    scenario,
    nirsevimab.avail = "yes"
  )
  # combine the datasets, the variable nirsevimab.avail denotes the two datasets
  rsv.data <- rbindlist(
    list(
      gen.rsv.hist,
      gen.rsv
    )
  )
  # convert data to factor, to simplify analysis
  rsv.data$adm.week <- factor(rsv.data$adm.week)
  rsv.data$nirsevimab <- factor(rsv.data$nirsevimab, levels = c("no", "yes"))
  rsv.data$nirsevimab.avail <- factor(rsv.data$nirsevimab.avail, levels = c("no", "yes"))
  rsv.data$age.grp <- factor(rsv.data$age.grp, levels = c("infant", "toddler", "child"))
  
  return(rsv.data)
}

do_analysis_prim <- function(birth.pop = 10000, vac.cov = 0, scenario = "maj") {
  # check we have correct vac.coverage. Efficacy analysis will not run with coverage 0 or 1
  # assertthat::assert_that(vac.cov > 0)
  # assertthat::assert_that(vac.cov < 1)
  
  # prepare containers to store results
  sig_prim <- c()
  sig_sec <- c()
  n_hosp <- c()
  vac_eff <- c()
  coef_nirse_eff <- c()
  
  # repeat analysis n.rep times
  for (i in 1:n.rep) {
    # sample vaccine efficacy from non-standard beta distribution using boundaries from published studies
    vac.eff <- extraDistr::rnsbeta(1, 6, 3, 0.5, 0.9)
    # prepare data with the current vac.eff
    data <- prep_analysis_data(birth.pop, vac.cov, vac.eff, scenario)
    # primary outcome = change in case numbers 
    model.prim <- glm(
      hosp ~ nirsevimab.avail + age.grp + adm.week,
      family = quasipoisson(link = "log"),
      data = data
    )
    # secondary outcome = vaccine efficacy estimation from hospital data
    # this only uses data where nirsevimab was available and compares children who received Nirsevimab to those that did not 
    # model may fail if there is no child who had nirsevimab or all children in hospital had nirsevimab. use try to avoid error
    try(model.vaceff <- glm(
      hosp ~ nirsevimab + age.grp + adm.week,
      family = quasipoisson(link = "log"),
      data = data[nirsevimab.avail == "yes", ]
    ), silent = TRUE)
    
    # store results
    sig_prim[i] <- tidy(model.prim)$p.value[2] <= sign.val
    n_hosp[i] <- data[nirsevimab == "yes", sum(hosp, na.rm = TRUE)]
    vac_eff[i] <- vac.eff
    sig_sec[i] <- ifelse(exists("model.vaceff"), tidy(model.vaceff)$p.value[2] <= sign.val, NA)
    coef_nirse_eff[i] <- ifelse(exists("model.vaceff"), 1 - exp(coef(model.vaceff)[2]), 0)
  }
  
  # return the data that we want from the models
  return(
    data.table(
      birth.pop = birth.pop,
      vac.cov = vac.cov,
      real.vac.eff = mean(vac_eff, na.rm = TRUE),
      real.vac.eff.l95ci = mean(vac_eff, na.rm = TRUE) - 1.96 * sd(vac_eff, na.rm = TRUE),
      real.vac.eff.u95ci = mean(vac_eff, na.rm = TRUE) + 1.96 * sd(vac_eff, na.rm = TRUE),
      est.vac.eff = mean(coef_nirse_eff, na.rm = TRUE),
      scenario = scenario,
      # av.hosp.rate = hosp.rate,
      power.prim = mean(sig_prim, na.rm = TRUE),
      power.sec = mean(sig_sec, na.rm = TRUE),
      n_hosp = round(mean(n_hosp, na.rm = TRUE))
    )
  )
}

# perform power analysis ####
pop <- c(10000, 15000, 20000, 25000, 30000)
vac.coverage <- seq(0, 1, 0.05)
scenario <- c("min", "maj")
params <- expand.grid(pop, vac.coverage, scenario)
names(params) <- c("cohort", "coverage", "scenario")

# this would likely benefit from running on parallel cores, could use the package furrr for this
container <- purrr::pmap(
  list(params$cohort, params$coverage, params$scenario),
  \(x, y, z) do_analysis_prim(birth.pop = x, vac.cov = y, scenario = z),
  .progress = list(name = "simulate RSV data", show_after = 10)  # this shows a progress bar, get your tea or coffee
)
container <- rbindlist(container)

# save the data ####
save(container, file = here("data", "rsv_insel_pwranalysis_nirsevimab_sim.RData"))

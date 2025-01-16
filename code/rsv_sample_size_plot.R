# Preloadare plots of power analysis and expectes Nirsevimab impact

library(here)
library(data.table)
library(ggplot2)
library(patchwork)
library(gt)

# load data ####
load(here("data", "rsv_insel_pwranalysis_nirsevimab_sim.RData"))

prep_plot_data <- function(
    birth.pop = 10000,
    vac.cov = 0,
    scenario = "maj"
) {
  plot.data <- data.table()
  # generate n repetitions of data
  for(i in 1:plot.rep) {
    vac.eff <- extraDistr::rnsbeta(1, 6, 3, 0.5, 0.9)
    rsv <- gen_rsv_data(
      birth.pop,
      vac.cov,
      vac.eff,
      scenario,
      nirsevimab.avail = "yes"
    )
    plot.data <- rbindlist(list(plot.data, rsv))
  }
  # calculate mean per week by dividing trough n repetitions
  plot.data <- plot.data[,
    .(n.hosp = round(sum(hosp) / plot.rep)),
    by = .(age.grp, adm.week)
  ]
  # add 0 for weeks without data
  plot.data <- plot.data[
    CJ(age.grp = unique(age.grp), adm.week = 1:52),
    on = c("age.grp", "adm.week")
  ]
  plot.data[, n.hosp := nafill(n.hosp, type = "const", fill = 0)]
  # add plot week to allow plotting of winter season
  plot.data[, year := fcase(
    adm.week < 27, plot.year + 1,
    default = plot.year
  )]
  # add date for pretty plotting
  plot.data[,
    date.season.plot := ISOweek::ISOweek2date(
      paste0(
        year, "-W",
        sprintf("%02d", adm.week),
        "-1"
      )
    )
  ]
  # change age.grp to factor
  plot.data$age.grp <- factor(plot.data$age.grp, levels = c("infant", "toddler", "child"))
  plot.data
}

# plots of power analysis - primary outcome ####
plots.prim <- list()
count <- 1
for (row in unique(container$birth.pop)) {
  for (col in unique(container$scenario)) {
    # define title of plot
    if(col == "maj") {
      title.txt <- paste0("Birth population: ", format(row, big.mark = "'"))
    }
    if(col == "min") {
      title.txt <- paste0("Birth population: ", format(row, big.mark = "'"))
    }
    p <- ggplot(container[birth.pop == row & scenario == col], aes(vac.cov, power.prim)) +
      geom_line(colour = "red", linewidth = 1.5) +
      geom_hline(yintercept = 0.9, linetype = 2) +
      scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), minor_breaks = NULL) +
      scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), minor_breaks = NULL) +
      labs(
        title = title.txt,
        x = "Immunisation coverage", y = "Power"
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey30"),
        text = element_text(size = 8)
      )
    plots.prim[[count]] <- p
    count <- count + 1
  }
}

p.prim <- wrap_elements(
  panel = grid::textGrob(
    paste0(
      "Scenario minor season\nAdmission rate infants: ",
      round(hosp.rate.min[age.grp == "infant", rate] * 100, 1), "%"
    ),
    gp = grid::gpar(fontsize = 10)
  )
) + 
  wrap_elements(
    grid::textGrob(
      paste0(
        "Scenario major season\nAdmission rate infants: ",
        round(hosp.rate.maj[age.grp == "infant", rate] * 100, 1), "%"
      ),
      gp = grid::gpar(fontsize = 10)
    )
  ) + 
  plots.prim[[1]] + plots.prim[[2]] + plots.prim[[3]] + plots.prim[[4]] + 
  plots.prim[[5]] + plots.prim[[6]] + plots.prim[[7]] + plots.prim[[8]] +
  plots.prim[[9]] + plots.prim[[10]] +
  plot_layout(ncol = 2, heights = c(1/11, rep(2/11, 5)))

# plots for power analysis - secondary outcome ####
plots.sec <- list()
count <- 1
for (row in unique(container$birth.pop)) {
  for (col in unique(container$scenario)) {
    # define title of plot
    if(col == "maj") {
      title.txt <- paste0("Birth population: ", format(row, big.mark = "'"))
    }
    if(col == "min") {
      title.txt <- paste0("Birth population: ", format(row, big.mark = "'"))
    }
    p <- ggplot(container[birth.pop == row & scenario == col & !is.na(power.sec)], aes(vac.cov, power.sec)) +
      geom_line(colour = "red", linewidth = 1.5) +
      geom_hline(yintercept = 0.9, linetype = 2) +
      scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), minor_breaks = NULL) +
      scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), minor_breaks = NULL) +
      labs(
        title = title.txt,
        x = "Immunisation coverage", y = "Power"
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey30"),
        text = element_text(size = 8)
      )
    plots.sec[[count]] <- p
    count <- count + 1
  }
}

p.sec <- wrap_elements(
  panel = grid::textGrob(
    paste0(
      "Scenario minor season\nAdmission rate infants: ",
      round(hosp.rate.min[age.grp == "infant", rate] * 100, 1), "%"
    ),
    gp = grid::gpar(fontsize = 10)
  )
) + 
  wrap_elements(
    grid::textGrob(
      paste0(
        "Scenario major season\nAdmission rate infants: ",
        round(hosp.rate.maj[age.grp == "infant", rate] * 100, 1), "%"
      ),
      gp = grid::gpar(fontsize = 10)
    )
  ) + plots.sec[[1]] + plots.sec[[2]] + plots.sec[[3]] + plots.sec[[4]] + 
  plots.sec[[5]] + plots.sec[[6]] + plots.sec[[7]] + plots.sec[[8]] +
  plots.sec[[9]] + plots.sec[[10]] + 
  plot_layout(ncol = 2, heights = c(1/11, rep(2/11, 5)))

# Plot bias we can expect for secondary analysis ####

p.bias.eff <- ggplot(
  container[birth.pop == 20000 & scenario == "maj" & vac.cov > 0 & vac.cov < 1, ],
  aes(vac.cov, est.vac.eff)
) +
  geom_ribbon(mapping = aes(ymin = real.vac.eff.l95ci, ymax = real.vac.eff.u95ci), fill = "grey70") +
  geom_line(mapping = aes(y = real.vac.eff)) +
  geom_line(colour = "red", linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), minor_breaks = NULL) +
  labs(
    x = "Immunisation coverage", y = "Estimated vaccine efficiency"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(colour = "grey90"),
    axis.line = element_line(colour = "grey30"),
    text = element_text(size = 8)
  )

# make plots for vaccination scenarios ####
make_plots <- function(birth.pop = 10000) {
  maj.season.nirse <- prep_plot_data(birth.pop, vac.cov = 0.5, scenario = "maj")
  maj.season.hist <- prep_plot_data(birth.pop, vac.cov = 0, scenario = "maj")
  min.season.nirse <- prep_plot_data(birth.pop, vac.cov = 0.5, scenario = "min")
  min.season.hist <- prep_plot_data(birth.pop, vac.cov = 0, scenario = "min")
  
  # make tables with overall numbers and percentages of age groups
  tab.maj.hist <- maj.season.hist[, .(n.hosp = sum(n.hosp)), by = age.grp] |> 
    _[, tot.hosp := sum(n.hosp)] |> 
    _[, perc := n.hosp / tot.hosp] |> 
    _[order(age.grp), ] |> 
    _[, .(age.grp, n.hosp, perc)] |> 
    _[, age.grp := c("< 1y", "1-2y", "> 2y")] |> 
    gt(rowname_col = "age.grp") |> 
    grand_summary_rows(columns = n.hosp, fns = list(All = ~ sum(., na.rm = TRUE))) |> 
    fmt_percent(columns = "perc", decimals = 0) |> 
    cols_label(age.grp = "Age", n.hosp = "N", perc = "%") |> 
    tab_options(table.font.size = 12)
    
  tab.maj.nirse <- maj.season.nirse[, .(n.hosp = sum(n.hosp)), by = age.grp] |> 
    _[, tot.hosp := sum(n.hosp)] |> 
    _[, perc := n.hosp / tot.hosp] |> 
    _[order(age.grp), ] |> 
    _[, .(age.grp, n.hosp, perc)] |> 
    _[, age.grp := c("< 1y", "1-2y", "> 2y")] |> 
    gt(rowname_col = "age.grp") |> 
    grand_summary_rows(columns = n.hosp, fns = list(All = ~ sum(., na.rm = TRUE))) |> 
    fmt_percent(columns = "perc", decimals = 0) |> 
    cols_label(age.grp = "Age", n.hosp = "N", perc = "%") |> 
    tab_options(table.font.size = 12)
  
  tab.min.hist <- min.season.hist[, .(n.hosp = sum(n.hosp)), by = age.grp] |> 
    _[, tot.hosp := sum(n.hosp)] |> 
    _[, perc := n.hosp / tot.hosp] |> 
    _[order(age.grp), ] |> 
    _[, .(age.grp, n.hosp, perc)] |> 
    _[, age.grp := c("< 1y", "1-2y", "> 2y")] |> 
    gt(rowname_col = "age.grp") |> 
    grand_summary_rows(columns = n.hosp, fns = list(All = ~ sum(., na.rm = TRUE))) |> 
    fmt_percent(columns = "perc", decimals = 0) |> 
    cols_label(age.grp = "Age", n.hosp = "N", perc = "%") |> 
    tab_options(table.font.size = 12)
  
  tab.min.nirse <- min.season.nirse[, .(n.hosp = sum(n.hosp)), by = age.grp] |> 
    _[, tot.hosp := sum(n.hosp)] |> 
    _[, perc := n.hosp / tot.hosp] |> 
    _[order(age.grp), ] |> 
    _[, .(age.grp, n.hosp, perc)] |> 
    _[, age.grp := c("< 1y", "1-2y", "> 2y")] |> 
    gt(rowname_col = "age.grp") |> 
    grand_summary_rows(columns = n.hosp, fns = list(All = ~ sum(., na.rm = TRUE))) |> 
    fmt_percent(columns = "perc", decimals = 0) |> 
    cols_label(age.grp = "Age", n.hosp = "N", perc = "%") |> 
    tab_options(table.font.size = 12)
  
  # make plots
  p.maj.hist <- ggplot(
    maj.season.hist[, .(n.hosp = sum(n.hosp)), by = .(date.season.plot)],
    aes(date.season.plot, n.hosp)
  ) +
    geom_col() +
    scale_x_date("", date_labels = "%b") +
    labs(title = "Major season RSV admissions\nwithout Nirsevimab", y = "n admitted to hospital") +
    theme(
      panel.background = element_rect(fill = "white")
    )
  p.maj.nirse <- ggplot(
    maj.season.nirse[, .(n.hosp = sum(n.hosp)), by = .(date.season.plot)],
    aes(date.season.plot, n.hosp)
  ) +
    geom_col() +
    scale_x_date("", date_labels = "%b") +
    labs(title = "Major season RSV admissions\nwith 50% Nirsevimab coverage", y = "n admitted to hospital") +
    theme(
      panel.background = element_rect(fill = "white")
    ) 
  
  p.min.hist <- ggplot(
    min.season.hist[, .(n.hosp = sum(n.hosp)), by = .(date.season.plot)],
    aes(date.season.plot, n.hosp)
  ) +
    geom_col() +
    scale_x_date("", date_labels = "%b") +
    labs(title = "Minor season RSV admissions\nwithout Nirsevimab", y = "n admitted to hospital") +
    theme(
      panel.background = element_rect(fill = "white")
    )
  p.min.nirse <- ggplot(
    min.season.nirse[, .(n.hosp = sum(n.hosp)), by = .(date.season.plot)],
    aes(date.season.plot, n.hosp)
  ) +
    geom_col() +
    scale_x_date("", date_labels = "%b") +
    labs(title = "Minor season RSV admissions\nwith 50% Nirsevimab coverage", y = "n admitted to hospital") +
    theme(
      panel.background = element_rect(fill = "white")
    ) 
  
  p.maj <- p.maj.hist +
    inset_element(tab.maj.hist, 0.1, 0.4, 0.5, 1) + 
    p.maj.nirse + 
    inset_element(tab.maj.nirse, 0.1, 0.4, 0.5, 1) &
    ylim(0, 30)
  
  p.min <- p.min.hist +
    inset_element(tab.min.hist, 0.1, 0.4, 0.5, 1) + 
    p.min.nirse + 
    inset_element(tab.min.nirse, 0.1, 0.4, 0.5, 1) &
    ylim(0, 30)
  
  list(p.maj, p.min)
}

# save plots ####
ggsave(
  here("figures", "rsvepi_power_age_prim.png"),
  plot = p.prim,
  scale = 0.75,
  width = 8, height = 12,
  dpi = 600,
)

ggsave(
  here("figures", "rsvepi_power_age_bias.png"),
  plot = p.bias.eff,
  scale = 0.6,
  width = 8, height = 6,
  dpi = 600
)

ggsave(
  here("figures", "rsvepi_power_age_sec.png"),
  plot = p.sec,
  scale = 0.75,
  width = 8, height = 12,
  dpi = 600,
)

# save figures ####
p.impact <- make_plots()
save(p.impact, file = here("data", "plots.impact.RData"))
save(list = c("container", "p.prim", "p.bias.eff", "p.sec"), file = here("data", "plots.pwr.RData"))

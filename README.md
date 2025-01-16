# rsv_sim
Repository contains the code to simulate RSV admission numbers at the Inselspital Bern.

Simulations are based on a historic dataset from the Inselspital Bern in Switzerland, ranging back to 1997.
Because of a relevant change in the number of children admitted per season, only data from summer 2008 onwards is used.
To be sure that the data reflect the typical biannual cycle of alternating minor and large season, only data until summer 2019 is included.
The dataset contains the following variables:

- adm.isow = week of hospital admission based on the ISO week date system.
To facilitate calculations there are always 52 weeks. Historical data from years with 53 ISO weeks was modified to include patients admitted in week 53 in week 52 and 1 of the following year (split 50/50 using 'sample()').
- age.grp = age at hospital admission categorised into three groups. Can be 'infant' (<= 365 days), 'toddler' (366-730 days) and 'child' (730-1826 days)
- season = strength of season. Can take the values 'min' or 'maj'
- season.nr = denotes an individual season (week 27 of previous year until week 26 of current year)

The code is written in R using the following packages:
'here' and 'data.table' as a backbone to do data wrangling.
'broom', 'extraDistr', and 'purrr' for the simulation.
'ggplot2', 'patchwork' and 'gt' to present results.

rm(list=ls(all=TRUE))
# data management sources
# part 1 - get forestry data into annual time-series
source("Data manage/Keogh forestry revised.R")
# part 2 - download climate data and integrate with other environmental data into annual time-series
source("Data manage/Get environmental data.R")
# part 3 - compile juvenile cohort that led to adult returns, adult cohort, and recruiting cohort assigned to brood year
source("Data manage/Data stock and recruit.R")
# part 4 - run MARSS DLM to recover missing environmental data from covariance of other variables
source("Data manage/Get missing data for environment.R")
# part 5 - run MARSS DLM to recover missing adult data from covariance with other species
source("Data manage/Get missing data for adults.R")
# part 6 - run MARSS DLM to recover missing outmigrating data from covariance with other species
source("Data manage/Get missing data for juv cohort.R")
# part 7 - run MARSS DLM to recover missing recruitment data from covariance with other species
source("Data manage/Get missing data for recruitment cohorts.R")
# part 8 - combine seal trend with NP salmon abundance trend to remove collinearity
source("Data manage/Combine collinear variables.R")
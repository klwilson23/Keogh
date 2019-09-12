source("some functions.R")
keogh <- readRDS("Keogh_stockRec_enviro.rds")
head(keogh)
pairs(keogh[,-grep(paste(c("pk","ct","co_","dv","ch_"),collapse="|",sep=""),colnames(keogh))])
columns <- c("sh_Smolts","sh_Adults","dv_Smolts","dv_Adults","ct_Smolts","ct_Adults","pk_Recruits","pk_Adults","co_Smolts","co_Adults","precip_3","win_precip_3","temp_3","npgo","mei","total","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

# steelhead
columns <- c("sh_Smolts","sh_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

# dolly varden
columns <- c("dv_Smolts","dv_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","total","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

# coastal cutthroat
columns <- c("ct_Smolts","ct_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

#coho
columns <- c("co_Smolts","co_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

#pink
columns <- c("pk_Recruits","pk_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

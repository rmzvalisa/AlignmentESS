

# Data and libraries

library(foreign)
library(car)
library(lavaan) 
library(plyr)
library(dplyr) 
library(readxl)
library(kableExtra)
library(MIE)
library(LittleHelpers)
library(gdata)
library(DataCombine)
library(ggplot2)
library(stringr)
library(tidyr)
library(gridExtra)
library(grid)
library(formattable)
source("ESS_functions.R")

# Downloading integrated and country specific data files, select and recode variables
## ESS 2 - country specific file for Italy
## ESS 3 - country specific file for Latvia and Romania
## ESS 4 - country specific file for Austria and Lithuania
## ESS 5 - country specific file for Austria 
## ESS 7 - country specific file for Russia
## ESS 9 - country specific file for Russia, Romania, Albania

## rlgatnd - How often attend religious services apart from special occasions
## pray - How often pray apart from at religious services
## rlgdgr - How religious are you

ESS_files <- c(
  "ESS1e06_6.sav", "ESS2e03_6.sav", "ESS2IT.sav", "ESS3e03_7.sav", "ESS3LV.sav",
  "ESS3RO.sav", "ESS4e04_5.sav", "ESS4LT.sav", "ESS4AT.sav", "ESS5e03_4.sav",
  "ESS5ATe1_1.sav", "ESS6e02_4.sav", "ESS7e02_2.sav", "Russian_social_survey_round_7.sav",
  "ESS8e02_1.sav", "ESS9e03.sav", "Russian_social_survey_round_9.sav", "ESS9ROe01.sav", 
  "ESS9ALe01.sav"
)

ESS_data <- lapply(ESS_files, read.spss, use.value.labels = F, to.data.frame = T)
for (i in 1:length(ESS_data)) {
  ESS_data[[i]] <- select(ESS_data[[i]], c(cntry, rlgatnd, pray, rlgdgr))
  colnames(ESS_data[[i]]) <- c("country", "attend", "pray", "person")
  
  # recode
  for (item in c("attend", "pray")) {
    ESS_data[[i]][, item] <- Recode(
      ESS_data[[i]][, item], rec = "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1; else=NA"
    )
  }
  ESS_data[[i]]$person = Recode(
    ESS_data[[i]]$person, rec = "0=1; 1=2; 2=3; 3=4; 4=5; 5=6; 6=7; 7=8; 8=9; 9=10; 10=11; else=NA"
  )
  
  ESS_data[[i]]$country <- as.factor(ESS_data[[i]]$country)
  
}

names(ESS_data) <- c(
  "ESS1", "ESS2_part", "ESS2_it", "ESS3_part", "ESS3_lv", "ESS3_ro",
  "ESS4_part", "ESS4_lt", "ESS4_at", "ESS5_part", "ESS5_at", "ESS6",
  "ESS7_part", "ESS7_ru", "ESS8", "ESS9_part", "ESS9_ru", "ESS9_ro", 
  "ESS9_al"
)

ESS_data$ESS2 <- rbind(ESS_data$ESS2_part, ESS_data$ESS2_it)
ESS_data$ESS3 <- rbind(ESS_data$ESS3_part, ESS_data$ESS3_lv, ESS_data$ESS3_ro)
ESS_data$ESS4 <- rbind(ESS_data$ESS4_part, ESS_data$ESS4_lt, ESS_data$ESS4_at)
ESS_data$ESS5 <- rbind(ESS_data$ESS5_part, ESS_data$ESS5_at)
ESS_data$ESS7 <- rbind(ESS_data$ESS7_part, ESS_data$ESS7_ru)
ESS_data$ESS9 <- rbind(ESS_data$ESS9_part, ESS_data$ESS9_ru, ESS_data$ESS9_ro, ESS_data$ESS9_al)

ESS_data <- ESS_data[names(ESS_data) %in%
                       c(
                         "ESS2_part", "ESS2_it", "ESS3_part", "ESS3_lv", "ESS3_ro", "ESS4_part",
                         "ESS4_lt", "ESS4_at", "ESS5_part", "ESS5_at", "ESS7_part", "ESS7_ru",
                         "ESS9_part", "ESS9_ru", "ESS9_ro", "ESS9_al"
                       ) ==
                       FALSE]

for (i in 1:length(ESS_data)) {
  ESS_data[[i]]$round <- names(ESS_data[i])
}

ESS_data <- ESS_data[order(names(ESS_data))]

# Clean
rm(ESS_files, item, i)


## share of missings
round(
  100 - 100/(nrow(do.call("rbind", ESS_data))/nrow(na.omit(do.call("rbind", ESS_data)))),
  0)

#-----------------------------------------------------------------------------------------------

# MGCFA
## remove Turkey from ESS 2 due to the low factor loading of frequency of religious attendance
ESS_data$ESS2 <- subset(ESS_data$ESS2, subset = !(ESS_data$ESS2$country %in% c("TR")))
ESS_data$ESS2$country <- droplevels(ESS_data$ESS2$country)


model_inv <- "Involv =~ person + attend + pray"

models_list <- lapply(
  ESS_data, function(x) globalMI.modified(
    model_inv, data = x, group = "country", estimator = "MLR", missing = "listwise"
  )
)
models_table <- rbind.fill(models_list)

colnames(models_table) <- c(
  "Model", "CFI", "Δ CFI", "RMSEA", "Δ RMSEA", "SRMR", "Δ SRMR", "CHI.sq",
  "DF"
)

models_table$CHI.sq <- as.numeric(as.character(models_table$CHI.sq)) %>%
  round(., 2)


kable(models_table) %>%
  group_rows("ESS 1", 1, 3) %>%
  group_rows("ESS 2", 4, 6) %>%
  group_rows("ESS 3", 7, 9) %>%
  group_rows("ESS 4", 10, 12) %>%
  group_rows("ESS 5", 13, 15) %>%
  group_rows("ESS 6", 16, 18) %>%
  group_rows("ESS 7", 19, 21) %>%
  group_rows("ESS 8", 22, 24) %>%
  group_rows("ESS 9", 25, 27) %>%
  footnote(
  general = "CFI = comparative fit index, RMSEA = root mean square error of approximation, 
           SRMR = standardized root mean residual, χ2 = chi-square, df = degrees of freedom.
           ESS 1 = European Social Survey Round 1, ESS 2 = European Social Survey Round 2, 
           ESS 3 = European Social Survey Round 3, ESS 4 = European Social Survey Round 4, 
           ESS 5 = European Social Survey Round 5, ESS 6 = European Social Survey Round 6, 
           ESS 7 = European Social Survey Round 7, ESS 8 = European Social Survey Round 8, 
           ESS 9 = European Social Survey Round 9."
  )



#-----------------------------------------------------------------------------------------------


# Alignment

## abbreviations for country numbers in alignment tables
codes_list <- lapply(ESS_data, function(x) as.data.frame(levels(x[, 1])))
for (i in 1:length(codes_list)) {
  codes_list[[i]] <- cbind(
    1:nrow(codes_list[[i]]),
    codes_list[[i]]
  )
  names(codes_list[i][[1]])[1] <- c("number")
  names(codes_list[i][[1]])[2] <- names(codes_list[i])
}
codes_list <- lapply(codes_list, function(x) trim(x))

codes_list <- Reduce(
  function(x, y) merge(x, y, by = "number", all = TRUE),
  codes_list
)
codes_list <- rbind(codes_list[10:32, ], codes_list[1:9, ])


for (i in 1:length(ESS_data)) {
  
  # create folders
  directory <- paste0("add your path to the folder here/", names(ESS_data[i]))
  folder <- dir.create(directory)
  setwd(directory)
  
  # run alignment
  runAlignment.modified(
    "Involv BY attend pray person", 
    group = "country", 
    dat = ESS_data[[i]],
    estim = "mlr", 
    sim.samples = c(100, 500, 1000, 1500, 2000),
    sim.reps = 500, 
    Mplus_com = "/Applications/Mplus/mplus", 
    summaries = F
  )
  
  # alignment main results
  ESS_align <- extractAlignment("fixed.out", silent = T)
  mv(
    "ESS_align", paste0(
      names(ESS_data[i]),
      "_align"
    )
  )
  
  # standard errors
  ESS_se <- extractSE.ESS("fixed.out")
  mv(
    "ESS_se", paste0(
      names(ESS_data[i]),
      "_se"
    )
  )
  
  
  # alignment simulation
  ESS_sim <- extractAlignmentSim(
    c("sim100.out", "sim500.out", "sim1000.out", "sim1500.out", "sim2000.out"),
    silent = T
  )
  ESS_sim <- cbind(
    ESS_sim[[1]][[1]][[1]], ESS_sim[[2]][[1]][[1]], ESS_sim[[3]][[1]][[1]],
    ESS_sim[[4]][[1]][[1]], ESS_sim[[5]][[1]][[1]]
  ) %>%
    as.data.frame(as.numeric(as.character(.)))
  
  colnames(ESS_sim) <- c("N = 100", "N = 500", "N = 1000", "N = 1500", "N = 2000")
  mv(
    "ESS_sim", paste0(
      names(ESS_data[i]),
      "_sim"
    )
  )
  
}

## re-estimate alignment for ESS4 - correct the group with factor mean = 0 from Cyprus to Latvia
## it is the closest one

setwd("add your path to the folder here/ESS4")

ESS4_align <- extractAlignment("fixed.out", silent = T)

ESS4_sim <- extractAlignmentSim(
  c("sim100.out", "sim500.out", "sim1000.out", "sim1500.out", "sim2000.out"),
  silent = T
) %>% 
  cbind(
  ESS4_sim[[1]][[1]][[1]], ESS4_sim[[2]][[1]][[1]], ESS4_sim[[3]][[1]][[1]],
  ESS4_sim[[4]][[1]][[1]], ESS4_sim[[5]][[1]][[1]]
) %>%
  as.data.frame(as.numeric(as.character(.)))
colnames(ESS4_sim) <- c("N = 100", "N = 500", "N = 1000", "N = 1500", "N = 2000")

ESS4_se <- extractSE.ESS("fixed.out")



ESS_align <- list(
  ESS1_align, ESS2_align, ESS3_align, ESS4_align, ESS5_align, ESS6_align,
  ESS7_align, ESS8_align, ESS9_align
)

# Figure with factor means across rounds
ESS_means <- ESS_align

ESS_se <- list(
  ESS1_se, ESS2_se, ESS3_se, ESS4_se, ESS5_se, ESS6_se, ESS7_se, ESS8_se,
  ESS9_se
)


round_names <- c("ESS 1", "ESS 2", "ESS 3", "ESS 4", "ESS 5", "ESS 6", "ESS 7", "ESS 8", "ESS 9")


for (i in 1:length(ESS_means)) {
  ESS_means[[i]] <- ESS_means[[i]][[4]][[1]][, 3:4] %>%
    na.omit()
  
  names(ESS_means[[i]])[names(ESS_means[[i]]) ==
                          "Factor.mean"] <- round_names[i]
  names(ESS_means[[i]])[names(ESS_means[[i]]) ==
                          "Group.value"] <- "code"
  
  code_ESS <- select(codes_list, c(number, names(ESS_data[i])))
  ESS_means[[i]] <- merge(
    ESS_means[[i]], code_ESS, by.x = c("code"),
    by.y = c("number"),
    all.x = TRUE
  )
  
  names(ESS_means[[i]])[2] <- round_names[i]
  names(ESS_means[[i]])[3] <- "country"
  
  ESS_means[[i]] <- cbind(ESS_means[[i]], ESS_se[[i]])
  ESS_means[[i]] <- ESS_means[[i]][, -1]
  ESS_means[[i]][, 1] <- formattable(ESS_means[[i]][, 1], digits = 2, format = "f")
  ESS_means[[i]][, 3] <- as.numeric(as.character(ESS_means[[i]][, 3]))
  names(ESS_means[[i]])[names(ESS_means[[i]]) ==
                          "ESS_se[[i]]"] <- "se"
}

plot1 <- ggmeans(ESS_means[[1]])
plot2 <- ggmeans(ESS_means[[2]])
plot3 <- ggmeans(ESS_means[[3]])
plot4 <- ggmeans(ESS_means[[4]])
plot5 <- ggmeans(ESS_means[[5]])
plot6 <- ggmeans(ESS_means[[6]])
plot7 <- ggmeans(ESS_means[[7]])
plot8 <- ggmeans(ESS_means[[8]])
plot9 <- ggmeans(ESS_means[[9]])


# pdf("Figure2_1.pdf",
#     width = 18.75, height = 27.08, 
#     bg = "transparent",      
#     colormodel = "cmyk")
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 2, ncol = 3)
# dev.off() 

# pdf("Figure2_2.pdf",
#     width = 18.75, height = 13.54, 
#     bg = "transparent",      
#     colormodel = "cmyk")
grid.arrange(plot7, plot8, plot9, nrow = 1, ncol = 3, 
             bottom = textGrob("          AL = Albania, AT = Austria, BE = Belgium, BG = Bulgaria, HR = Croatia, CY = Cyprus, CZ = Czech Republic, DK = Denmark, EE = Estonia, FI = Finland,
            FR = France, DE = Germany, GR = Greece, HU = Hungary, IS = Iceland, IE = Ireland, IL = Israel, IT = Italy, XK = Kosovo, LV = Latvia,  LT = Lithuania, 
           LU = Luxembourg, ME = Montenegro, NL = Netherlands, NO = Norway, PL = Poland, PT = Portugal, RO = Romania, RU = Russian Federation, RS = Serbia, 
       SK = Slovakia, SI = Slovenia, ES = Spain, SE = Sweden, CH = Switzerland, TR = Turkey, UA = Ukraine, GB = United Kingdom",
                               gp = gpar(fontface = 3, fontsize = 17.5)))
# dev.off() 





# Alignment tables

## extract specific fit indices
for (i in 1:length(ESS_align)) {
  ESS_align[[i]] <- ESS_align[[i]][[3]][, c(3, 6, 4, 5)]
  ESS_align[[i]] <- cbind(
    c(
      "attend.int", "person.int", "pray.int", "attend.load", "person.load", "pray.load"
    ),
    ESS_align[[i]]
  )
  colnames(ESS_align[[i]]) <- c(
    "", "R-Square", "Fit Contribution", "Invariant Groups", "Non-invariant Groups"
  )
  name <- names(ESS_data[i])
  ESS_align[[i]] <- FindReplace(
    data = ESS_align[[i]], Var = "Invariant Groups", replaceData = codes_list,
    from = "number", to = as.character(name),
    exact = F, vector = F
  )
  ESS_align[[i]] <- FindReplace(
    data = ESS_align[[i]], Var = "Non-invariant Groups", replaceData = codes_list,
    from = "number", to = as.character(name),
    exact = F, vector = F
  )
}



# Share of non-invariant countries

share_ni <- ESS_align
for (i in 1:length(share_ni)) {
  ni_countries <- sapply(share_ni[[i]][5], function(x) strsplit(x, " ")) %>%
    lapply(., length) %>%
    do.call("rbind", .)
  
  total_countries <- sapply(share_ni[[i]][4], function(x) strsplit(x, " ")) %>%
    lapply(., length) %>%
    do.call("rbind", .)
  total_countries <- sum(ni_countries[1] + total_countries[1])
  
  share_ni[[i]] <- cbind(share_ni[[i]][1], ni_countries)
  share_ni[[i]][2] <- round(share_ni[[i]][2]/total_countries, 2)
  
  names(share_ni[[i]])[names(share_ni[[i]]) ==
                             "Var.1"] <- "parameter"
  names(share_ni[[i]])[names(share_ni[[i]]) ==
                             "ni_countries"] <- paste0("NI share_", names(ESS_data[i]))
  
}

share_ni <- Reduce(
  function(x, y) merge(x, y, by = "parameter", all = TRUE),
  share_ni
)

# R-square

rsquare <- ESS_align
for (i in 1:length(rsquare)) {
  rsquare[[i]] <- rsquare[[i]][, 1:2]
  names(rsquare[[i]])[1] <- "parameter"
  rsquare[[i]][, 2] <- as.numeric(as.character(rsquare[[i]][, 2])) %>%
    round(., 2)
  
  names(rsquare[[i]])[names(rsquare[[i]]) ==
                        "R-Square"] <- paste0("R-Square_", names(ESS_data[i]))
}

rsquare <- Reduce(
  function(x, y) merge(x, y, by = "parameter", all = TRUE),
  rsquare
)


# Fit function contribution

fit_function <- ESS_align
for (i in 1:length(fit_function)) {
  fit_function[[i]] <- fit_function[[i]][, c(1, 3)]
  names(fit_function[[i]])[1] <- "parameter"
  names(fit_function[[i]])[2] <- paste0("Fit-Function_", names(ESS_data[i]))
  fit_function[[i]][, 2] <- as.numeric(as.character(fit_function[[i]][, 2])) %>%
    round(., 2)
}

fit_function <- Reduce(
  function(x, y) merge(x, y, by = "parameter", all = TRUE),
  fit_function
)


# Share of countries by zones

Country_abb <- read_excel("Country_abb.xlsx")

share_zones <- ESS_align
for (i in 1:length(share_zones)) {
  share_zones[[i]] <- FindReplace(
    data = share_zones[[i]], Var = "Invariant Groups", replaceData = Country_abb,
    from = "abbreviation", to = "zone", exact = F, vector = F
  )
  share_zones[[i]] <- FindReplace(
    data = share_zones[[i]], Var = "Non-invariant Groups", replaceData = Country_abb,
    from = "abbreviation", to = "zone", exact = F, vector = F
  )
  zones <- sapply(share_zones[[i]][4], function(x) strsplit(x, " "))
  Catholic_inv <- str_count(zones, "Catholic")
  Protestant_inv <- str_count(zones, "Protestant")
  Orthodox_inv <- str_count(zones, "Orthodox")
  Islamic_inv <- str_count(zones, "Islamic")
  Israel_inv <- str_count(zones, "Israel")
  
  zones_total <- strsplit(share_zones[[i]][1, 5], " ")
  Catholic_total <- Catholic_inv[1] + str_count(zones_total, "Catholic")[1]
  Protestant_total <- Protestant_inv[1] + str_count(zones_total, "Protestant")[1]
  Orthodox_total <- Orthodox_inv[1] + str_count(zones_total, "Orthodox")[1]
  Islamic_total <- Islamic_inv[1] + str_count(zones_total, "Islamic")[1]
  Israel_total <- Israel_inv[1] + str_count(zones_total, "Israel")[1]
  
  
  Catholic_inv <- round(Catholic_inv/Catholic_total, 2)
  Protestant_inv <- round(Protestant_inv/Protestant_total, 2)
  Orthodox_inv <- round(Orthodox_inv/Orthodox_total, 2)
  Islamic_inv <- round(Islamic_inv/Islamic_total, 2)
  Israel_inv <- round(Israel_inv/Israel_total, 2)
  
  share_zones[[i]] <- cbind(
    share_zones[[i]][1], Catholic_inv, Protestant_inv, Orthodox_inv,
    Islamic_inv, Israel_inv
  )
  colnames(share_zones[[i]]) <- c(
    "parameter", paste0("Catholic_", names(ESS_data[i])),
    paste0("Protestant_", names(ESS_data[i])),
    paste0("Orthodox_", names(ESS_data[i])),
    paste0("Islamic_", names(ESS_data[i])),
    paste0("Israel_", names(ESS_data[i]))
  )
  
  
}

## Ignore the warnings, they are related to the arguments' format only

share_zones <- Reduce(
  function(x, y) merge(x, y, by = "parameter", all = TRUE),
  share_zones
)


## TABLE A2 IN APPENDIX

# Combine into one table

align_table <- Reduce(
  function(x, y) merge(x, y, all = TRUE),
  list(rsquare, fit_function, share_ni, share_zones)
)
align_table <- align_table %>%
  gather(variable, value, -parameter) %>%
  mutate(
    location = sub(
      ".*(ESS1|ESS2|ESS3|ESS4|ESS5|ESS6|ESS7|ESS8|ESS9).*", "\\1",
      variable
    ),
    variable = sub("_?(ESS1|ESS2|ESS3|ESS4|ESS5|ESS6|ESS7|ESS8|ESS9)_?", "", variable)
  ) %>%
  spread(variable, value)

align_table <- cbind(
  align_table[, 1:2], align_table[, 10], align_table[, 4], align_table[, 7], 
  align_table[, 3], align_table[, 8:9], align_table[, 5:6]
)

align_table$parameter[duplicated(align_table$parameter)] <- NA
align_table[, 2:10][is.na(align_table[, 2:10])] <- "--"

colnames(align_table) <- c(
  "Parameter", "Round", "R-square", "Fit function", "NI %", "Catholic",
  "Orthodox", "Protestant", "Islamic", "Israel"
)

align_table$`Fit function` <- round(align_table$`Fit function`, 1)

kable(align_table) %>%
  add_header_above(
    c(
      ` ` = 1, ` ` = 1, `Fit indices` = 2, ` ` = 1, `Share of invariant countries, by zone` = 5
    )
  ) %>%
  group_rows("Frequency of religious attendance", 1, 18) %>%
  group_rows("Self-assessed religiosity", 19, 36) %>%
  group_rows("Frequency of praying", 37, 54) %>%
  footnote(
    general = "Fit Function = Fit Function Contribution, 
  NI % = Share of non-invariant parameters, 
  -- = countries did not participate in survey. 
  For the ESS rounds abbreviations, see the note to Table 1.
  Baseline group and N groups: 
           ESS 1 - Czech Republic, 22; ESS 2 - Czech Republic, 25; 
           ESS 3 - Germany, 25; ESS 4 - Latvia, 31; ESS 5 - Estonia, 28; 
           ESS 6 - Czech Republic, 29; ESS 7 -  Czech Republic, 22; 
           ESS 8 -  Czech Republic, 23; ESS 9 -  Czech Republic, 32"
  )



# Main table

## mean for share of NI countries
share_ni$Share <- rowMeans(share_ni[, 2:10])
## round(colMeans(share_ni[, 2:10]), 2)

## mean for R-square
rsquare$Rsquare <- rowMeans(rsquare[, 2:10])
## round(colMeans(rsquare[, 2:10]), 2)

align_table_sum <- as.data.frame(cbind(share_ni$Share, rsquare$Rsquare))
align_table_sum <- round(align_table_sum, 2)


## mean for Fit function contribution

fit_function[nrow(fit_function) + 1, ] <- NA
fit_function[7, 2:10] <- colMeans(fit_function[1:6, 2:10])
fit_function[7, 1] <- "Mean"

## replace the differences for the fit function contribution with the model specific mean
fit_function[1:6, 2] <- fit_function[7, 2] - fit_function[1:6, 2]
fit_function[1:6, 3] <- fit_function[7, 3] - fit_function[1:6, 3]
fit_function[1:6, 4] <- fit_function[7, 4] - fit_function[1:6, 4]
fit_function[1:6, 5] <- fit_function[7, 5] - fit_function[1:6, 5]
fit_function[1:6, 6] <- fit_function[7, 6] - fit_function[1:6, 6]
fit_function[1:6, 7] <- fit_function[7, 7] - fit_function[1:6, 7]
fit_function[1:6, 8] <- fit_function[7, 8] - fit_function[1:6, 8]
fit_function[1:6, 9] <- fit_function[7, 9] - fit_function[1:6, 9]
fit_function[1:6, 10] <- fit_function[7, 10] - fit_function[1:6, 10]

fit_function[, 2:10] <- round(fit_function[, 2:10], 1)


align_table_sum <- as.data.frame(
  cbind(fit_function[1:6, 1], align_table_sum, fit_function[1:6, 2:10])
)

colnames(align_table_sum) <- c("", "", "", round_names)

kable(align_table_sum) %>%
  add_header_above(
    c(
      `Parameter` = 1, `NI %` = 1, `R-square` = 1, `Fit Function Contribution` = 9
    )
  ) %>%
  group_rows("Frequency of religious attendance", 1, 2) %>%
  group_rows("Self-assessed religiosity", 3, 4) %>%
  group_rows("Frequency of praying", 5, 6) %>%
  footnote(
    general = "NI % = Share of non-invariant parameters.
    R2 and the share of noninvariant countries are the averages across nine survey rounds.
    The fit contribution is presented separately for each survey because it is not standardized; 
    thus, its averaging across surveys would be biased. 
    The positive values indicate the contribution higher than its average value in a survey, 
    while the negative values indicate the contribution lower than its average value in a survey. 
    The means at the bottom of the table are the raw means for each survey.
  For the ESS rounds abbreviations, see the note to Table 1.
  Baseline group and N groups: 
           ESS 1 - Czech Republic, 22; ESS 2 - Czech Republic, 25; 
           ESS 3 - Germany, 25; ESS 4 - Latvia, 31; ESS 5 - Estonia, 28; 
           ESS 6 - Czech Republic, 29; ESS 7 -  Czech Republic, 22; 
           ESS 8 -  Czech Republic, 23; ESS 9 -  Czech Republic, 32"
  )





# Simulations

ESS_sim <- rbind(
  ESS1_sim, ESS2_sim, ESS3_sim, ESS4_sim, ESS5_sim, ESS6_sim, ESS7_sim,
  ESS8_sim, ESS9_sim
)

ESS_sim <- apply(ESS_sim, 2, function(x)
  round(
    as.numeric(as.character(x)), 3
  ))

rownames(ESS_sim) <- names(ESS_data)
kable(ESS_sim) %>%
  add_header_above(c(Round = 1, `Number of observations per group` = 5)) %>%
  footnote(
    general = "For the ESS rounds abbreviations, see the note to Table 1."
  )



## TABLE A1 IN APPENDIX

obs_total <- lapply(ESS_data, function(x) na.omit(x))
obs_total <- lapply(obs_total, function(x) trim(x))

for (i in 1:length(obs_total)) {
  obs_total[[i]] <- data.frame(table(obs_total[[i]]$country))
  names(obs_total[i][[1]])[1] <- c("Country")
  names(obs_total[i][[1]])[2] <- names(obs_total[i])
}

obs_total <- Reduce(
  function(x, y) merge(x, y, by = "Country", all = TRUE),
  obs_total
)

obs_total <- trim(obs_total)
obs_total <- obs_total[order(obs_total[, 1]), ]

obs_total <- rbind(
  obs_total, data.frame(Country = "Total", t(colSums(obs_total[, -1], na.rm = T)))
)
obs_total <- rbind(
  obs_total, data.frame(Country = "N countries", t(colSums(!is.na(obs_total[-39, -1]))))
)
obs_total[is.na(obs_total)] <- "--"


obs_total <- merge(obs_total, Country_abb[, 2:3], by.x = "Country", by.y = "abbreviation", all.x = T)
kable(obs_total) %>%
  footnote(
  general = "AL = Albania, AT = Austria, BE = Belgium, BG = Bulgaria, HR = Croatia, 
           CY = Cyprus, CZ = Czech Republic, DK = Denmark, EE = Estonia, FI = Finland, 
           FR = France, DE = Germany, GR = Greece, HU = Hungary, IS = Iceland, IE = Ireland, 
           IL = Israel, IT = Italy, XK = Kosovo, LV = Latvia,  LT = Lithuania, LU = Luxembourg,
           ME = Montenegro, NL = Netherlands, NO = Norway, PL = Poland, PT = Portugal, RO = Romania,
           RU = Russian Federation, RS = Serbia, SK = Slovakia, SI = Slovenia, ES = Spain, 
           SE = Sweden, CH = Switzerland, TR = Turkey, UA = Ukraine, GB = United Kingdom.
           Catholic = Catholic Europe zone, Islamic = African-Islamic zone, 
           Orthodox = Orthodox Europe zone, Protestant = Protestant Europe zone.
           ESS 1 = European Social Survey Round 1, ESS 2 = European Social Survey Round 2, 
           ESS 3 = European Social Survey Round 3, ESS 4 = European Social Survey Round 4, 
           ESS 5 = European Social Survey Round 5, ESS 6 = European Social Survey Round 6, 
           ESS 7 = European Social Survey Round 7, ESS 8 = European Social Survey Round 8, 
           ESS 9 = European Social Survey Round 9.
           Turkey was excluded from the ESS 2 sample due to the low factor loading of 
           the frequency of religious attendance indicator."
  )






# Session info ----------------------------------------------------

sessionInfo(package = NULL)

# attached packages: formattable_0.2.1 gridExtra_2.3 tidyr_1.1.3 stringr_1.4.0 DataCombine_0.2.21  
# gdata_2.18.0 LittleHelpers_0.5-10 lme4_1.1-27.1 Matrix_1.3-2 magrittr_2.0.1      
# MIE_0.5-3 shinyWidgets_0.6.0 shinyjs_2.0.0 igraph_1.2.6 DT_0.17             
# shiny_1.6.0 ggforce_0.3.3 ggrepel_0.9.1 ggplot2_3.3.5 kableExtra_1.3.4    
# readxl_1.3.1 dplyr_1.0.7 plyr_1.8.6 lavaan_0.6-7 car_3.0-10          
# carData_3.0-4 oreign_0.8-81      



library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(minpack.lm)
library(broom)
library(stringr)

tresp_data <- read_csv("curve_level_data.csv")
head(tresp_data)


#the quadratic term

tresp_data <- tresp_data%>%
  mutate(Tleaf_mean2 = Tleaf_mean^2)
head(tresp_data)

#create new plant ids
tresp_data <- tresp_data %>%
  mutate(plant_id_new = str_extract(id, "^[^_]+"))
head(tresp_data)

#fitting model per plant id for vcmax

quad_exp_fit_lm <- tresp_data %>%
  group_by(plant_id_new) %>% #group by plant id
  filter(!is.na(obs_vcmax)) %>%  # remove NA rows
  filter(n() >= 3) %>%                               # require at least 3 points
  group_modify(~{
    model <- lm(log(obs_vcmax) ~ Tleaf_mean + Tleaf_mean2, data = .x)
    data.frame(
      n_points = nrow(.x),
      intercept = coef(model)[1],
      b1 = coef(model)[2],
      b2 = coef(model)[3],
      Topt = -coef(model)[2] / (2 * coef(model)[3])
    )
  })
head(quad_exp_fit_lm)
summary(quad_exp_fit_lm)

#writing the model fit as excel sheet
library(writexl)

write_xlsx(quad_exp_fit_lm, "vcmax_tempresp_fits.xlsx")


#fitting model per plant id for jcmax

quad_jmax_fit_lm <- tresp_data %>%
  group_by(plant_id_new) %>% #group by plant id
  filter(!is.na(obs_jmax)) %>%  # remove NA rows
  filter(n() >= 3) %>%                               # require at least 3 points
  group_modify(~{
    model <- lm(log(obs_jmax) ~ Tleaf_mean + Tleaf_mean2, data = .x)
    data.frame(
      n_points = nrow(.x),
      intercept = coef(model)[1],
      b1 = coef(model)[2],
      b2 = coef(model)[3],
      Topt = -coef(model)[2] / (2 * coef(model)[3])
    )
  })
head(quad_jmax_fit_lm)

#writing the model fit as excel sheet

write_xlsx(quad_jmax_fit_lm, "jmax_tempresp_fits.xlsx")

#upload extracted a, b, and c values with treatment information for vcmax 
vcmax_tempresp_fits <- read_csv("vcmax_tempresp_fits.csv")
head(vcmax_tempresp_fits)
view(vcmax_tempresp_fits)

#check distributions for a, b, and c
hist(vcmax_tempresp_fits$ a, main = "Distribution of a")
hist(vcmax_tempresp_fits$ b, main = "Distribution of b")
hist(vcmax_tempresp_fits$ c, main = "Distribution of c")

#using only stable fits
vcmax_tempresp_fits <- vcmax_tempresp_fits %>%
  filter(c < 0,b <= 10, a <=10)
head(vcmax_tempresp_fits)

#check distribution for good fits
hist(vcmax_tempresp_fits$ a, main = "Distribution of a")
hist(vcmax_tempresp_fits$ b, main = "Distribution of b")
hist(vcmax_tempresp_fits$ c, main = "Distribution of c")



##fit lme model for parameter a
a_lmer <- lmer(a ~ Nfert * airtemp_factor + 
                  (1|Rack:airtemp_factor), data = vcmax_tempresp_fits)
plot(residuals(a_lmer)~fitted(a_lmer))
summary(a_lmer)
Anova(a_lmer)
emmeans(a_lmer, ~ Nfert * airtemp_factor)

##fit lme model for parameter b
b_lmer <- lmer(b ~ Nfert * airtemp_factor + 
                 (1|Rack:airtemp_factor), data = vcmax_tempresp_fits)
plot(residuals(b_lmer)~fitted(b_lmer))
summary(b_lmer)
Anova(b_lmer)
emmeans(b_lmer, ~ Nfert * airtemp_factor)

##fit lme model for parameter c
c_lmer <- lmer(c ~ Nfert * airtemp_factor + 
                 (1|Rack:airtemp_factor), data = vcmax_tempresp_fits)
plot(residuals(c_lmer)~fitted(c_lmer))
summary(c_lmer)
Anova(c_lmer)
emmeans(c_lmer, ~ Nfert * airtemp_factor)

#writing the good fits
write_xlsx(vcmax_tempresp_fits, "good_fits.xlsx")

#upload extracted a, b, and c values with treatment information for jmax
jmax_tempresp_fits <- read_csv("jmax_tempresp_fits.csv")
head(jmax_tempresp_fits)

#chlorophyll calculations

#upload chlorophyll data
data <- read_csv("chlorophyll extraction.csv")
head(data)
view(data)
data$avg_649 <- (data$A649 + data$A649.1 + data$A649.2)/3
data$avg_665 <- (data$A665 + data$A665.1 + data$A665.2)/3
data$chla_ug.ml <- (12.19 * data$avg_665) - (3.45 * data$avg_649) # ug mL-1, from Wellburn (1994)
data$chlb_ug.ml <- (21.99 * data$avg_649) - (5.32 * data$avg_665) # ug mL-1, from Wellburn (1994)
data$chla_g.ml <- data$chla_ug.ml / 1000000 # g mL-1
data$chlb_g.ml <- data$chlb_ug.ml / 1000000 # g mL-1
data$chla_g <- data$chla_g.ml * 10 # 10 mL of DMSO
data$chlb_g <- data$chlb_g.ml * 10 # 10 mL of DMSO
data$chla_g.m2 <- data$chla_g / (data$chl_area / 10000) # convert area to m2
data$chlb_g.m2 <- data$chlb_g / (data$chl_area / 10000) # convert area to m2
data$chla_mol.m2 <- data$chla_g.m2 / 893.51 # 893.51 g mol-1 chlorophyll a
data$chlb_mol.m2 <- data$chlb_g.m2 / 907.47 # 907.47 g mol-1 chlorophyll b
data$chla_mmol.m2 <- data$chla_mol.m2 * 1000
data$chlb_mmol.m2 <- data$chlb_mol.m2 * 1000





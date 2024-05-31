# Install necessary packages if not present
# install.packages(c("broom", "tidyr", "dplyr", "purrr", "readr", "here"))

# Load packages
library(broom)
library(tidyr)
library(dplyr)
library(purrr)
library(readr)
library(readxl)
library(stringr)
library(tidyverse)

#Load file
Sphingolipids_Lpg2_to_CW <- read_excel("C:/Users/micfr287/Work Folders/Desktop/Sphingolipids Lpg2 to CW.xlsx")
View(Sphingolipids_Lpg2_to_CW)
Sphingonew <- Sphingolipids_Lpg2_to_CW
View(Sphingonew)

#Set Palmitate as a reference group
Sphingonew$Group <- as.factor(Sphingonew$Group)
Sphingonew <- mutate(Sphingonew,
                     diet = ifelse(Group == "Palmitate", 0,
                                   ifelse(Group == "Linoleate", 1, "NA")))

#Rename sex
Sphingonew <- Sphingonew %>%
  rename(sex = `sex (1=kvinna, 2=man)`)

#Rename all lipids and apoB
Sphingonew <- Sphingonew %>%
  rename(apob_pre = apoB_w0)
Sphingonew <- Sphingonew %>%
  rename(apob_post = apoB_w8)
Sphingonew <- Sphingonew %>%
  rename(sm_fa14_0_pre = SM140_w0)
Sphingonew <- Sphingonew %>%
  rename(sm_fa14_0_post = SM140_w8)
Sphingonew <- Sphingonew %>%
  rename(cer_fa18_0_pre = Cer180_w0)
Sphingonew <- Sphingonew %>%
  rename(cer_fa18_0_post = Cer180_w8)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa18_0_pre = DiCer180_w0)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa18_0_post = DiCer180_w8)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa20_0_pre = DiCer200_w0)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa20_0_post = DiCer200_w8)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa22_0_pre = DiCer220_w0)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa22_0_post = DiCer220_w8)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa24_1_pre = DiCer241_w0)
Sphingonew <- Sphingonew %>%
  rename(dicer_fa24_1_post = DiCer241_w8)
Sphingonew <- Sphingonew %>%
  rename(glucer_fa18_0_pre = GluCer180_w0)
Sphingonew <- Sphingonew %>%
  rename(glucer_fa18_0_post = GluCer180_w8)

# Ceramides expressed in micromole/L instead of nanomole/L

Sphingonew$cer_fa18_0_pre <- Sphingonew$cer_fa18_0_pre/1000
Sphingonew$cer_fa18_0_post <- Sphingonew$cer_fa18_0_post/1000
Sphingonew$dicer_fa18_0_pre <- Sphingonew$dicer_fa18_0_pre/1000
Sphingonew$dicer_fa18_0_post <- Sphingonew$dicer_fa18_0_post/1000
Sphingonew$dicer_fa20_0_pre <- Sphingonew$dicer_fa20_0_pre/1000
Sphingonew$dicer_fa20_0_post <- Sphingonew$dicer_fa20_0_post/1000
Sphingonew$dicer_fa22_0_pre <- Sphingonew$dicer_fa22_0_pre/1000
Sphingonew$dicer_fa22_0_post <- Sphingonew$dicer_fa22_0_post/1000
Sphingonew$dicer_fa24_1_pre <- Sphingonew$dicer_fa24_1_pre/1000
Sphingonew$dicer_fa24_1_post <- Sphingonew$dicer_fa24_1_post/1000
Sphingonew$glucer_fa18_0_pre <- Sphingonew$glucer_fa18_0_pre/1000
Sphingonew$glucer_fa18_0_post <- Sphingonew$glucer_fa18_0_post/1000
Sphingonew$sm_fa14_0_pre <- Sphingonew$sm_fa14_0_pre/1000
Sphingonew$sm_fa14_0_post <- Sphingonew$sm_fa14_0_post/1000

# Define a function to calculate a weighted sphingolipid score based on specific sphingolipids
# The weights are derived from the observed intervention effect in DIVAS


# scale score by observed total effect in DIVAS
make_sl_score <- function(.sm_fa14_0, .cer_fa18_0, .dicer_fa18_0, .dicer_fa20_0, .dicer_fa22_0, .dicer_fa24_1, .glucer_fa18_0) {
  score <- (log(.sm_fa14_0)*-0.16 +
              log(.cer_fa18_0)*-0.168 +
              log(.dicer_fa18_0)*-0.148 +
              log(.dicer_fa20_0)*-0.133 +
              log(.dicer_fa22_0)*-0.125 +
              log(.dicer_fa24_1)*-0.122 +
              log(.glucer_fa18_0)*-0.17)/0.355
  
  return(score)
}


# z-scaled score
make_sl_score_z <- function(.sm_fa14_0, .cer_fa18_0, .dicer_fa18_0, .dicer_fa20_0, .dicer_fa22_0, .dicer_fa24_1, .glucer_fa18_0) {
  score <- scale((log(.sm_fa14_0)*-0.16 +
                    log(.cer_fa18_0)*-0.168 +
                    log(.dicer_fa18_0)*-0.148 +
                    log(.dicer_fa20_0)*-0.133 +
                    log(.dicer_fa22_0)*-0.125 +
                    log(.dicer_fa24_1)*-0.122 +
                    log(.glucer_fa18_0)*-0.17)) |>
    as.numeric()
  
  return(score)
}


### Prerequisites for Analysis:

# 1. Handle Missing Values: Remove participants with missing data or use standard imputation methods.
# 2. Adjust Variable Names: Ensure lipid and other variable names (age, sex, BMI, APOB) in the script match your dataset. Modify as needed.
# 3. Adapt Input Dataset: Ensure the script points to the correct input dataset file or database.
# 4. Adapt Output Path: Modify the script to specify the desired output path for saving results.


# Please adapt the following variable names :
# diet
# bmi
# age
# sex
# sphingolipids
# pre_ and post_intervention_indicator string


pre_intervention_indicator <- "_pre"
post_intervention_indicator <- "_post"


lipid_names <- Sphingonew |>
  dplyr::select(contains("cer") |
                  dplyr::contains("sm")) |>
  names() |>
  stringr::str_remove(paste0(pre_intervention_indicator,
                             "|",
                             post_intervention_indicator)) |>
  unique()


lm.wrapper <- function(.df, .var) {
  # rename variables to unify
  temp <- .df %>%
    dplyr::rename(pre = str_c(!!.var, pre_intervention_indicator),
                  post = str_c(!!.var, post_intervention_indicator)) |>
    dplyr::mutate(pre = log(pre),
                  post = log(post))
  
  model <- lm(post ~ pre + diet + pre + age + bmi + sex,
              data = temp)
  
  return(model)
}

lm.wrapper_z <- function(.df, .var) {
  # rename variables to unify
  temp <- .df %>%
    dplyr::rename(pre = str_c(!!.var, pre_intervention_indicator),
                  post = str_c(!!.var, post_intervention_indicator)) |>
    dplyr::mutate(pre = log(pre) |> scale() |> as.numeric(),
                  post = log(post)|> scale() |> as.numeric())
  
  model <- lm(post ~ pre + diet + pre + age + bmi + sex,
              data = temp)
  
  return(model)
}


lipid_results <- tibble::tibble(lipid = c(lipid_names, "apob")) |>
  dplyr::mutate(res = purrr::map(lipid, ~lm.wrapper(Sphingonew, .x)),
                res_z = purrr::map(lipid, ~lm.wrapper_z(Sphingonew, .x)),
                coeffs = purrr::map(res, broom::tidy),
                coeffs_z = purrr::map(res_z, broom::tidy),
                model_stats = purrr::map(res, broom::glance),
                model_stats_z = purrr::map(res_z, broom::glance),
                residuals = purrr::map(res, resid),
                residuals_z = purrr::map(res_z, resid)) |>
  dplyr::select(-c(res, res_z)) # need to remove model objects, because they contain original data


# Unnest results to check results
lipid_results |> tidyr::unnest(coeffs)
lipid_results |>  tidyr::unnest(model_stats)

# Define the output path and save the final_results object
output_path <- here::here("C:/Users/micfr287/Work Folders/Desktop/Andra forskningsprojekt utanför PhD/lipid_results.rds") # adapt output path
readr::write_rds(x = lipid_results, file = output_path)


input_df <- Sphingonew|>
  # Calculate sphingolipid scores before and after the intervention using the defined function
  dplyr::mutate(score_pre = make_sl_score(.cer_fa18_0  = cer_fa18_0_pre,
                                          .dicer_fa18_0 = dicer_fa18_0_pre,
                                          .dicer_fa20_0 = dicer_fa20_0_pre,
                                          .dicer_fa22_0 = dicer_fa22_0_pre,
                                          .dicer_fa24_1 = dicer_fa24_1_pre,
                                          .glucer_fa18_0 = glucer_fa18_0_pre,
                                          .sm_fa14_0 = sm_fa14_0_pre),
                score_post = make_sl_score(.cer_fa18_0  = cer_fa18_0_post,
                                           .dicer_fa18_0 = dicer_fa18_0_post,
                                           .dicer_fa20_0 = dicer_fa20_0_post,
                                           .dicer_fa22_0 = dicer_fa22_0_post,
                                           .dicer_fa24_1 = dicer_fa24_1_post,
                                           .glucer_fa18_0 = glucer_fa18_0_post,
                                           .sm_fa14_0 = sm_fa14_0_post),
                score_pre_z = make_sl_score_z(.cer_fa18_0  = cer_fa18_0_pre,
                                              .dicer_fa18_0 = dicer_fa18_0_pre,
                                              .dicer_fa20_0 = dicer_fa20_0_pre,
                                              .dicer_fa22_0 = dicer_fa22_0_pre,
                                              .dicer_fa24_1 = dicer_fa24_1_pre,
                                              .glucer_fa18_0 = glucer_fa18_0_pre,
                                              .sm_fa14_0 = sm_fa14_0_pre),
                score_post_z = make_sl_score_z(.cer_fa18_0  = cer_fa18_0_post,
                                               .dicer_fa18_0 = dicer_fa18_0_post,
                                               .dicer_fa20_0 = dicer_fa20_0_post,
                                               .dicer_fa22_0 = dicer_fa22_0_post,
                                               .dicer_fa24_1 = dicer_fa24_1_post,
                                               .glucer_fa18_0 = glucer_fa18_0_post,
                                               .sm_fa14_0 = sm_fa14_0_post))


# Perform linear modeling to assess the intervention effect, including demographic and baseline variables
score_results <- tibble::tibble(model = c("sl_1", "sl_2", "sl_1_z", "sl_2_z"), # Model identifiers
                                # Define models with and without baseline and post-intervention APOB levels
                                res = list(lm(score_post ~ diet + score_pre + bmi + age + sex, data = input_df),
                                           lm(score_post ~ diet + score_pre + bmi + age + sex + apob_pre + apob_post, data = input_df),
                                           lm(score_post_z ~ diet + score_pre_z + bmi + age + sex, data = input_df),
                                           lm(score_post_z ~ diet + score_pre_z + bmi + age + sex + apob_pre + apob_post, data = input_df)),
                                # Extract model statistics and coefficients
                                model_stats = purrr::map(res, broom::glance),
                                coeffs = purrr::map(res, broom::tidy),
                                residuals = purrr::map(res, resid)) |>
  dplyr::select(-res) # need to remove model objects, because they contain original data


# Unnest results to check results
score_results |> tidyr::unnest(coeffs)
score_results |>  tidyr::unnest(model_stats)


# Define the output path and save the final_results object
output_path <- here::here("C:/Users/micfr287/Work Folders/Desktop/Andra forskningsprojekt utanför PhD/score_results.rds") # adapt output path
readr::write_rds(x = score_results, file = output_path)


# Missing rows
ids_diet_missing <- Sphingonew |>
  select(diet, apob_pre, sm_fa14_0_post) |>
  mutate(across(c(apob_pre, sm_fa14_0_post), is.na), # indicate which values are missing, which are present. I assume it will be the same participant across all sphingolipids
         id = row_number())

ids_diet_missing

# Define the output path and save the final_results object 
output_path <- here::here("C:/Users/micfr287/Work Folders/Desktop/Andra forskningsprojekt utanför PhD/info.rds") # adapt output path 
readr::write_rds(x = ids_diet_missing, file = output_path)


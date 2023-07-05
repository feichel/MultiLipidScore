


#' give informative names to food groups
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
rename_food_groups <- function(df) {
  df %>%
  rename("whole grain bread" = group1,
         "other bread" = group2,
         "grain flakes, grains, muesli" = group3,
         "cornflakes, crisps" = group4,
         "pasta, rice" = group5 ,
         "vegetarian dishes" = group6 ,
         "chips" = group7 ,
         "pizza" = group8 ,
         "cake, cookies" = group9 ,
         "confectionary" = group10,
         "sweet bread spreads" = group11,
         "eggs" = group12,
         "fresh fruit" = group13,
         "canned fruit" = group14,
         "raw vegetables" = group15,
         "cabbage" = group16,
         "cooked vegetables" = group17,
         "garlic" = group18,
         "mushrooms" = group19,
         "legumes" = group20,
         "potatoes" = group21,
         "fried potatoes" = group22,
         "nuts" = group23,
         "low-fat dairy products" = group24,
         "high-fat dairy products" = group25,
         "low-fat cheese" = group26,
         "high-fat cheese" = group27,
         "water" = group28,
         "coffee" = group29,
         "de-caf coffee" = group30,
         "tea" = group31,
         "fruit juice" = group32,
         "low-energy soft drinks" = group33,
         "high-energy soft drinks" = group34,
         "beer" = group35,
         "wine" = group36,
         "spirits" = group37,
         "other alcoholic beverages" = group38,
         "butter" = group39,
         "margarine" = group40,
         "other vegetable fat" = group41,
         "other fat" = group42,
         "sauce" = group43,
         "deserts" = group44,
         "fish" = group45,
         "poultry" = group46,
         "meat" = group47,
         "processed meat" = group48,
         "soup" = group49) %>%
    janitor::clean_names()
}




#' Rename lipids classes (i.e. "base") in accordance with lipidmaps
#'
#' @param .class_var
#'
#' @return
#' @export
#'
#' @examples
reorder_base_new <- function(.class_var) {
  order <- c("CE",
             "FFA",
             "MG",
             "DG",
             "TG",
             "SM",
             "Cer",
             "dhCer",
             "LacCer",
             "HexCer",
             "PI" ,
             "PC",
             "LPC",
             "PE",
             "PEP",
             "PEO",
             "LPE")

  forcats::fct_relevel(.class_var, order)
}



#' Function to run linear regression to assess effect of diet on baseline-adjusted
#' post intervention differences between diets
#'
#' @param .df
#' @param .var
#'
#' @return
#' @export
#'
#' @examples
lm.wrapper <- function(.df, .var) {
  # rename variables to unify
  temp <- .df %>%
    rename(pre = str_c(!!.var, "_v1"),
           post = str_c(!!.var, "_v2"))

  model <- lm(post ~ diet + pre + age + bmi + sex,
              data = temp)

  return(model)
}



#' Log-transform and scale continuous variables. SD and mean are taken from sub-cohort
#' and are applied to sub-cohort memebers and external cases
#'
#' @param .df
#' @param .var
#'
#' @return
#' @export
#'
#' @examples
scale_subcohort <- function(.df, .var) {

  var_filtered <- .df %>%
    filter(ia_subcohort == 1) %>%
    pull({{.var}})

  mean <- var_filtered %>% log %>% mean

  sd <- var_filtered %>% log %>% sd

  var <- .df %>%
    pull({{.var}})

  return((log(var) - mean)/sd)
}








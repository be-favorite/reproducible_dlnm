# 0. loading libraries --------------------------------------------------
library(tidyverse)
library(fpp3)
library(forecast)
library(dlnm)
library(splines)
library(foreach)
library(parallel)
library(doParallel)
library(furrr)
library(mvmeta)
library(latex2exp)
library(scales)
library(rlist) # list.save(), list.load()
library(patchwork)
ggplot2::theme_set(theme_classic())

# 1. ready to covariates name --------------------------------------------
covariates <- chicagoNMMAPS |> 
    as_tibble() |> 
    select(temp, rhum, dptp)

# 2. preparing custom functions for optimizing DLNMs ----------------------
## (1) making formula
outcome <- "n"
exposure <- "cb"
weekly_season <- paste(c("s1_7", "c1_7"), collapse = " + ")
daily_season <- paste(c("s1_365", "c1_365", "s2_365", "c2_365", "s3_365", "c3_365"), 
                      collapse = " + ")
get_formula <- function(.covariates){
    if(length(.covariates) == 0){
        return(paste(outcome,
                     paste(exposure, daily_season, weekly_season, sep = " + "),
                     sep = " ~ ") %>% 
                   as.formula)
    }else if(length(.covariates) == 1){
        return(paste(outcome,
                     paste(exposure, daily_season, weekly_season, .covariates, sep = " + "),
                     sep = " ~ ") %>% 
                   as.formula)
    }else{
        return(paste(outcome,
                     paste(exposure, daily_season, weekly_season, paste(.covariates, collapse = " + "), 
                           sep = " + "),
                     sep = " ~ ") %>% 
                   as.formula)
    }
}
combine.lists <- function(list1, list2){
    # Combine lists 'list1' and 'list2', giving precedence to elements found in 'list2':
    # that is, if $something is found in both 'list1' and 'list2',
    # the new (output) list will have the same values as 'list2' in $something
    
    # Version 1.0 (August 2017)
    #
    # Function developed by
    # Patrick Belisle
    # Division of Clinical Epidemiology
    # McGill University Hospital Center
    # Montreal, Qc, Can
    #
    # patrick.belisle@rimuhc.ca
    # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
    
    
    list1.names <- names(list1)
    list2.names <- names(list2)
    
    new.list <- list1
    
    
    tmp <- match(list2.names, list1.names)
    w <- which(!is.na(tmp))
    
    if (length(w) > 0)
    {
        # take values from list2 in matching dimension names
        tmp <- tmp[!is.na(tmp)]
        new.list[[tmp]] <- list2[[w]]
        
        # append elements of 'list2' with unmatched names
        new.list <- c(new.list, list2[-w])
    }
    else
    {
        new.list <- c(new.list, list2)
    }
    
    new.list
}
covariates_comb <- combine.lists( # Model lists for the best subset selection
    combn(colnames(covariates)[1:3], 0, simplify = FALSE), 
    combn(colnames(covariates)[1:3], 1, simplify = FALSE)
) %>% 
    combine.lists(
        combn(colnames(covariates)[1:3], 2, simplify = FALSE)
    ) %>% 
    combine.lists(
        combn(colnames(covariates)[1:3], 3, simplify = FALSE)
    )
names(covariates_comb) <- 1:(2^ncol(covariates)) # 2^k -> k is the number of covariates you consider

### (2) for Tuning what covariates include, maximum lag, df of var dimension, df of lag dimension
fqaic <- function(model) { # Q-AIC FUNCTION
    loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
    phi <- summary(model)$dispersion
    qaic <- -2*loglik + 2*summary(model)$df[3]*phi
    return(qaic)
}
dlnm_model <- function(.data, .pollution, .index, .formula){ # function for modeling
    cb <<- crossbasis(.pollution, 
                      lag = grid[[.index, ".lag"]], 
                      argvar = list(fun = "ns", df = grid[[.index, "var_df"]]),
                      arglag = list(fun = "ns", df = grid[[.index, "lag_df"]]))
    m <- glm(formula = .formula, family = quasipoisson(), data = .data)
    fqaic(m)
}
grid <- expand.grid(.lag = 7:31, var_df = 2:5, lag_df = 2:5) %>% # Grid for hyper-parameters
    as_tibble()
dlnm_tuned <- function(.d, .p, .f){
    grid %>% 
        mutate(qaic = future_map_dbl(seq(nrow(grid)), 
                                     ~dlnm_model(.data = .d, 
                                                 .pollution = .p, 
                                                 .index = .x,
                                                 .formula = .f))) %>% 
        arrange(qaic)
}
dlnm_newyork <- function(.tb, airpollution, .by, .at){ # First, optimizing a model for Newyork
    mydata <- .tb %>% filter(city == "Newyork") # get data
    
    # do parallel processing
    numCores <- parallel::detectCores() - 1
    myCluster <- parallel::makeCluster(numCores)
    ## register it to be used by %dopar%
    doParallel::registerDoParallel(cl = myCluster)
    tuned <- foreach::foreach(i = seq_along(covariates_comb), .combine = rbind) %dopar% {
        source("./paper_reproducible_tutorial_custom_function.R")
        future::plan(cluster, workers = parallel::detectCores() - 1)
        dlnm_tuned(mydata, mydata %>% pull(airpollution), get_formula(covariates_comb[[i]])) %>% 
            slice(1) %>% 
            as.numeric()
    }
    # Finally, it is always recommendable to stop the cluster when we are done working with it.
    parallel::stopCluster(cl = myCluster)
    
    tuned_best <- tuned %>%
        as_tibble() %>% 
        rename(.lag = V1, var_df = V2, lag_df = V3, qaic = V4) %>% 
        mutate(formula_number = 1:8) %>%
        arrange(qaic, formula_number)
    tuned_best2 <- tuned_best %>% 
        mutate(diff = qaic - tuned_best$qaic[1]) %>% 
        # the algorithm selects the simplest model 
        # if there are models with a QAIC diffrences of less than 2 from the optimal model.
        filter(diff < 2)
    
    .lag <- tuned_best2 %>% 
        pull(.lag) %>% 
        min()
    var_df <- tuned_best2 %>% 
        pull(var_df) %>% 
        min()
    lag_df <- tuned_best2 %>% 
        pull(lag_df) %>% 
        min()
    formula_number <- tuned_best2 %>% 
        pull(formula_number) %>% 
        min()
    f <- get_formula(covariates_comb[[formula_number]])
    cb <<- crossbasis(mydata %>% pull(airpollution), # get cross-basis function for the single pollutant model
                      lag = .lag,
                      argvar = list(df = var_df),
                      arglag = list(df = lag_df))
    .model <- glm(f, family = quasipoisson(), data = mydata) # fitting the best models
    pred <- crosspred(cb, .model,
                      cen = .tb %>% filter(city == "Newyork") %>% pull(airpollution) |> median(), 
                      by = .by, at = .at)
    
    highlow_effect <- .tb %>% 
        filter(city == "Newyork") %>% 
        pull(airpollution) %>% 
        quantile(c(0.1, 0.9)) %>% 
        round(1)
    pred2 <- crosspred(cb, .model,
                       cen = .tb %>% filter(city == "Newyork") %>% pull(airpollution) |> median(), 
                       by = .by, at = c(highlow_effect[1], highlow_effect [2]))
    
    final_list <- list(tuned_best, f, cb, .model, pred, pred2) # 모형에 관한 정보 저장
    names(final_list) <- c("tuned_best", "formula", "cb", "model_spec", "predictions", "highlow")
    list.save(final_list, 
              str_c("./Best models/", airpollution, "_dlnm_Newyork", ".RData"))
}
### All hyper-parameters in the models except Newyork were selected based on the model for Newyork
### to do multivariate meta-analysis.
dlnm_model2 <<- function(.data, .pollution, .formula, .model){ # function for modeling
    cb <<- crossbasis(.pollution, 
                      lag = attr(.model$cb, which  = "lag")[2], # maximum lag days
                      argvar = list(fun = "ns", df = attr(.model$cb, which  = "df")[1]),
                      arglag = list(fun = "ns", df = attr(.model$cb, which  = "df")[2]))
    m <- glm(formula = .formula, family = quasipoisson(), data = .data)
    fqaic(m)
}
dlnm_best <- function(.tb, .city, airpollution, .model, .by, .at){
    mydata <- .tb %>% filter(city == .city)
    future::plan(cluster, workers = parallel::detectCores() - 1)
    
    tuned_best <- covariates_comb %>% 
        map(~get_formula(.x)) %>% 
        future_map(~dlnm_model2(.data = mydata,
                                .pollution = mydata %>% pull(airpollution),
                                .formula = .x,
                                .model = .model)) %>% 
        as_tibble() %>% 
        pivot_longer(1:8, names_to = "formula_number", values_to = "qaic") %>% 
        mutate(formula_number = as.numeric(formula_number)) %>% 
        arrange(qaic, formula_number)
    formula_number <- tuned_best %>% 
        mutate(diff = qaic - tuned_best$qaic[1]) %>% 
        filter(diff < 2) %>% 
        pull(formula_number) %>% 
        min()
    f <- get_formula(covariates_comb[[formula_number]])
    cb <<- crossbasis(mydata %>% pull(airpollution), # get cross-basis function for the single pollutant model
                      lag = attr(.model$cb, which  = "lag")[2], 
                      argvar = list(fun = "ns", df = attr(.model$cb, which  = "df")[1]),
                      arglag = list(fun = "ns", df = attr(.model$cb, which  = "df")[2]))
    .model <- glm(f, family = quasipoisson(), data = mydata) # fitting the best models
    highlow_effect <- .tb %>% 
        filter(city == .city) %>% 
        pull(airpollution) %>% 
        quantile(c(0.1, 0.9)) %>% 
        round(1)
    
    pred <- crosspred(cb, .model,
                      # Newyork 기준 centering
                      cen = .tb %>% filter(city == "Newyork") %>% pull(airpollution) |> median(), 
                      by = .by, at = .at)
    pred2 <- crosspred(cb, .model,
                       # Newyork 기준 centering
                       cen = .tb %>% filter(city == "Newyork") %>% pull(airpollution) |> median(), 
                       by = .by, at = c(highlow_effect[1], highlow_effect [2]))
    
    final_list <- list(tuned_best, f, cb, .model, pred, pred2)
    names(final_list) <- c("tuned_best", "formula", "cb", "model_spec", "predictions", "highlow")
    list.save(final_list, 
              str_c("./Best models/", airpollution, "_dlnm_", .city, ".RData"))
}
get_model <- function(.city, airpollution){
    list.load(str_c("./Best models/", airpollution, "_dlnm_", .city, ".RData"))
}
# 3. preparing multivariate meta-analysis --------------------------------
meta <- function(.tb, airpollution, .model, .by, .at){
    ####################################################################
    # Preparing and fitting Multivariate meta-analysis
    ####################################################################
    
    # Study area
    city <- .tb |>  
        pull(city) |> 
        unique()
    
    m <- length(city)
    df <- prod(attr(.model$cb, which  = "df"))
    
    ymat <- matrix(NA, nrow = m, ncol = df,
                   dimnames = list(city, str_c("spl", seq(df))))
    Slist <-  vector("list", m)
    names(Slist) <- city
    
    for (i in seq(m)){
        .model <- get_model(city[[i]], airpollution)
        
        cat(i, "")
        ymat[i,] <- .model$predictions$coef
        Slist[[i]] <- .model$predictions$vcov
    }
    
    # MULTIVARIATE META-ANALYSIS (fixed-effect)
    mv <- mvmeta(ymat, Slist, method="ml")
    
    ####################################################################
    # CREATE BASIS FOR PREDICTION
    ####################################################################
    
    # A cross basis matrix for the multivariate meta-analysis is created with
    # the PM10 average of 3 cities.
    air_avg <- .tb %>% 
        group_by(date) %>% 
        summarize(
            pm10 = mean(pm10),
        )
    cb <- crossbasis(air_avg %>% pull(airpollution),
                     lag = attr(.model$cb, which  = "lag")[2],
                     argvar = list(df = attr(.model$cb, which  = "df")[1]),
                     arglag = list(df = attr(.model$cb, which  = "df")[2]))
    cen <- air_avg %>%
        select(-date) %>% 
        map_dfr(~median(.x)) %>% 
        pull(airpollution)
    highlow_effect <- air_avg %>% 
        pull(airpollution) %>% 
        quantile(c(0.1, 0.9)) %>% 
        round(1)
    
    ####################################################################
    # PREDICTION FROM MODELS
    ####################################################################
    
    list("predictions" = crosspred(cb, coef = coef(mv), vcov = vcov(mv), model.link="log",
                                   by = .by, at = .at, 
                                   cen = cen),
         "highlow" = crosspred(cb, coef = coef(mv), vcov = vcov(mv), model.link="log",
                               by = .by, at = c(highlow_effect[1], highlow_effect[2]), 
                               cen = cen))
    
}

# 4.  Preparing for visualization -----------------------------------------
plot_3d <- function(.city, airpollution){
    .model <- get_model(.city, airpollution)
    .model[["predictions"]] %>% 
        plot(xlab = toupper(airpollution), zlab = "RR", 
             main = .city, 
             theta = 210, phi = 30, lphi = 30,
             border = "gray40", ticktype = "detailed", font.main = 1, 
             cex.main = 2, cex.lab = 1.5, cex.axis = 1.2, 
             cex.sub = 1.5)
    
}

plot_3d_meta <- function(.meta, airpollution){
    .meta[["predictions"]] %>% 
        plot(xlab = toupper(airpollution), zlab = "RR", 
             theta = 210, phi = 30, lphi = 30,
             border = "gray40", ticktype = "detailed", font.main = 1, 
             cex.main = 2, cex.lab = 1.5, cex.axis = 1.2, 
             cex.sub = 1.5)
}

plot_overall <- function(.tb, .city, airpollution, .xlab){
    .model <- get_model(.city, airpollution)
    .model[["predictions"]] %>% 
        plot(col = "tomato", lwd = 2, "overall",
             ylab = "RR", xlab = .xlab, 
             main = .city, 
             font.main = 1, cex.main = 2.5, cex.lab = 2.2, cex.axis = 2,
             cex.sub = 1.5,
             sub = str_c("Maximum lag: ", attr(.model$cb, "lag")[2], "days"))
    .tb %>% 
        filter(city == .city) %>% 
        pull(airpollution) %>% 
        rug(., quiet = TRUE)
}

plot_overall_meta <- function(.meta, .tb, airpollution, .xlab){
    .meta[["predictions"]] %>% 
        plot(col = "tomato", lwd = 2, "overall",
             ylab = "RR", xlab = .xlab, 
             font.main = 1, cex.main = 2.5, cex.lab = 2.2, cex.axis = 2,
             cex.sub = 1.5,
             sub = str_c("Maximum lag: ", .meta[["predictions"]]$lag[2], "days"))
    .tb %>% 
        group_by(date) %>% 
        summarize(
            pm10 = mean(pm10)
        ) %>% 
        pull(airpollution) %>% 
        rug(., quiet = TRUE)
}

make_png <- function(airpollution, w, h, type){
    png(str_c("./plot/", "reproducible_", airpollution, "_", type, ".png"),
        res = 100, width = w, height = h)
}

plot_slice <- function(.tb, .city, airpollution, math_air){
    .model <- get_model(.city, airpollution)
    pred <- .model[["highlow"]] 
    highlow_effect <- .tb %>% 
        filter(city == .city) %>% 
        pull(airpollution) %>% 
        quantile(c(0.1, 0.9)) %>% 
        round(1) %>% 
        as.character()
    
    effect_labs <- c(TeX(str_c("Low-", "$", math_air ,"$", " effect")),
                     TeX(str_c("High-", "$", math_air ,"$", " effect")))
    tibble(
        lag = rep(seq(pred$lag[[1]], pred$lag[[2]]), 2),
        effect = rep(c("Low", "High"), 
                     each = length(seq(pred$lag[[1]], pred$lag[[2]]))),
        RR = c(pred$matRRfit[highlow_effect[[1]], ], pred$matRRfit[highlow_effect[[2]], ]),
        low_RR = c(pred$matRRlow[highlow_effect[[1]], ], pred$matRRlow[highlow_effect[[2]], ]),
        high_RR = c(pred$matRRhigh[highlow_effect[[1]], ], pred$matRRhigh[highlow_effect[[2]], ])
    ) %>% 
        mutate(.col = ifelse(low_RR > 1 & high_RR > 1, "+",
                             ifelse(low_RR < 1 & high_RR < 1, "-", "0")),
               effect = factor(effect, levels = c("Low", "High"), labels = effect_labs)) %>% 
        ggplot(., aes(x = lag, y = RR)) +
        geom_errorbar(aes(x = lag, y = RR, ymin = low_RR, ymax = high_RR)) +
        geom_point(aes(col = .col), shape = 19, size = 2) +
        geom_hline(yintercept = 1) +
        scale_y_continuous(breaks= pretty_breaks()) +
        scale_x_continuous(breaks= pretty_breaks()) +
        scale_color_manual(values = c("+" = "red", "0" = "grey40", "-" = "blue")) +
        labs(x = "Lag",
             title = .city) + 
        theme(text = element_text(size = 20), legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~ effect, labeller = "label_parsed")
}

plot_slice_meta <- function(.tb, .meta, airpollution, math_air){
    pred <- .meta[["highlow"]]
    highlow_effect <- .tb %>% 
        group_by(date) %>% 
        summarize(
            pm10 = mean(pm10)
        ) %>% 
        pull(airpollution) %>% 
        quantile(c(0.1, 0.9)) %>% 
        round(1) %>% 
        as.character()
    
    effect_labs <- c(TeX(str_c("Low-", "$", math_air ,"$", " effect")),
                     TeX(str_c("High-", "$", math_air ,"$", " effect")))
    tibble(
        lag = rep(seq(pred$lag[[1]], pred$lag[[2]]), 2),
        effect = rep(c("Low", "High"), 
                     each = length(seq(pred$lag[[1]], pred$lag[[2]]))),
        RR = c(pred$matRRfit[highlow_effect[[1]], ], pred$matRRfit[highlow_effect[[2]], ]),
        low_RR = c(pred$matRRlow[highlow_effect[[1]], ], pred$matRRlow[highlow_effect[[2]], ]),
        high_RR = c(pred$matRRhigh[highlow_effect[[1]], ], pred$matRRhigh[highlow_effect[[2]], ])
    ) %>% 
        mutate(.col = ifelse(low_RR > 1 & high_RR > 1, "+",
                             ifelse(low_RR < 1 & high_RR < 1, "-", "0")),
               effect = factor(effect, levels = c("Low", "High"), labels = effect_labs)) %>% 
        ggplot(., aes(x = lag, y = RR)) +
        geom_errorbar(aes(x = lag, y = RR, ymin = low_RR, ymax = high_RR)) +
        geom_point(aes(col = .col), shape = 19, size = 2) +
        geom_hline(yintercept = 1) +
        scale_y_continuous(breaks= pretty_breaks()) +
        scale_x_continuous(breaks= pretty_breaks()) +
        scale_color_manual(values = c("+" = "red", "0" = "grey40", "-" = "blue")) +
        labs(x = "Lag") + 
        theme(text = element_text(size = 20), legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~ effect, labeller = "label_parsed")
}

make_png2 <- function(airpollution, w, h){
    ggsave(str_c("./plot/", "reproducible_", airpollution ,"_highlow", ".png"),
           width = w, height = h, dpi = 250)
}
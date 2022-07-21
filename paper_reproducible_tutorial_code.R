# 1. Loading custom function --------------------------------------------------
source("./paper_reproducible_tutorial_custom_function.R")
ggplot2::theme_set(theme_classic())

# 2. Load and being ready to data -----------------------------------------------------------
chicago <- chicagoNMMAPS |> 
    as_tibble() |> 
    select(date, resp, pm10, temp, rhum, dptp) |> 
    mutate(city = "Chicago",
           pm10 = ifelse(is.na(pm10) | pm10 < 0 | pm10 > 150, mean(pm10, na.rm = TRUE), pm10),
           rhum = ifelse(is.na(rhum), mean(rhum, na.rm = TRUE), rhum)) |> 
    rename(N = resp)
# outcome: resp(Respiratory deaths in chicago from 1987 - 2000)
# exposure: pm10
# seasonality: fourier series (we're going to add in third step)
# candidate of covariates: temp, rhum, dptp
newyork <- chicago |> 
    mutate(city = "Newyork")
miami <- chicago |> 
    mutate(city = "Miami")

adjust <- function(.data, add_n, m_pm, sd_pm, m_temp, sd_temp,
                   m_rhum, sd_rhum, m_dptp, sd_dptp){
    .data$N <- .data$N + rpois(n = nrow(.data), lambda = add_n)
    .data$pm10 <- .data$pm10 + rnorm(n = nrow(.data), mean = m_pm, sd = sd_pm)
    .data$temp <- .data$temp + rnorm(n = nrow(.data), mean = m_temp, sd = sd_temp)
    .data$rhum <- .data$rhum + rnorm(n = nrow(.data), mean = m_rhum, sd = sd_rhum)
    .data$dptp <- .data$dptp + rnorm(n = nrow(.data), mean = m_dptp, sd = sd_dptp)
    .data |> 
        mutate(
            rhum = ifelse(rhum > 100, 100, rhum)
        )
}
set.seed(1234)
chicago2 <- adjust(chicago, add_n = 10, m_pm = 15, sd_pm = 3, m_temp = 1, sd_temp = 1,
       m_rhum = 0.1, sd_rhum = 0.1, m_dptp = 1, sd_dptp = 1)
newyork2 <- adjust(newyork, add_n = 20, m_pm = 20, sd_pm = 5, m_temp = -1, sd_temp = 1.5,
                   m_rhum = 1, sd_rhum = 0.5, m_dptp = 2, sd_dptp = 0.1)
miami2 <- adjust(miami, add_n = 3, m_pm = 5, sd_pm = 1, m_temp = 3, sd_temp = 0.5,
                   m_rhum = 4, sd_rhum = 0.5, m_dptp = 3, sd_dptp = 0.5)

usa_death <- bind_rows(
    chicago2, newyork2, miami2
)

seasonality <- msts(usa_death %>% 
                        filter(city == "Chicago"), 
                    seasonal.periods = c(7, 365.25), 
                    start = c(1987, 1)) %>% 
    fourier(K = c(1, 3)) %>% 
    as_tibble %>% 
    mutate(date = seq(ymd("1987-01-01"), ymd("2000-12-31"), by = "day")) %>% 
    select(date, everything())

# 3. EDA -----------------------------------------------------------------
order_meteor <- c("temp", "dptp", "rhum")
labels_meteor <- c("temp" = "Temperature (C°)",
                   "dptp" = "Dew point temperature (C°)",
                   "rhum" = "Relative humidity (%rh)")

## Box-plots for airpollution
boxplot_air <- usa_death |> 
    pivot_longer(pm10) |> 
    ggplot(aes(x = city, y = value)) +
    geom_boxplot() +
    labs(x = "", y = TeX("$PM_{10} \\, (\\mu g /m^3)$")) +
    theme(text = element_text(family = "sans", size = 15),
          legend.position = "none") +
    coord_flip()

## Box-plots for meterological factors
boxplot_meteor <- usa_death |> 
    pivot_longer(temp:dptp) |> 
    mutate(
        name = factor(name, levels = order_meteor, labels = labels_meteor)
    ) |> 
    ggplot(aes(x = city, y = value)) +
    geom_boxplot() +
    labs(x = "", y = "") +
    facet_wrap(~ name, scales = "free", ncol = 1) +
    theme(text = element_text(family = "sans", size = 15),
          legend.position = "none") +
    coord_flip()

## monthly incidences plot
plot_monthly <- usa_death |> 
    mutate(date = yearmonth(date)) |> 
    group_by(date, city) |> 
    summarize(N = sum(N)) |> 
    ggplot(aes(x = date, y = N)) +
    geom_line() +
    facet_wrap(~ city, scales = "free", ncol = 1) +
    labs(x = "", y = "") +
    theme(text = element_text(family = "sans", size = 15),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none")

## multiple plots with {patchwork}
plot_monthly + (boxplot_air / boxplot_meteor) +
    plot_layout(width = c(1.5, 2)) +
    plot_annotation(tag_levels = "A")

# 4. Add fourier series for seasonalities ----------------------------
## add 2 fourier terms for weekly seasonalities, 
## and add 6 fourier terms for daily seasonalities.
seasonality <- msts(usa_death %>% 
                        filter(city == "Chicago"), 
                    seasonal.periods = c(7, 365.25), 
                    start = c(1987, 1)) %>% 
    fourier(K = c(1, 3)) %>% 
    as_tibble %>% 
    mutate(date = seq(ymd("1987-01-01"), ymd("2000-12-31"), by = "day")) %>% 
    select(date, everything())
usa_death2 <- usa_death %>% 
    left_join(seasonality, by = "date") %>%
    janitor::clean_names()

# 5. Optimizing DLNM for each cities and doing multivariate meta-analysis ---------------------------------------------------------
## (1) PM10
dlnm_newyork(usa_death2, "pm10", 10, 20:90)
.model <- get_model("Newyork", "pm10")
dlnm_best(usa_death2, "Miami", "pm10", .model, 10, 20:90)
dlnm_best(usa_death2, "Chicago", "pm10", .model, 10, 20:90)
meta_pm10 <- meta(usa_death2, "pm10", .model, 10, 20:90)

# 6. Visualization -------------------------------------------------------
make_png("pm10", 1800, 1000, "3D")
par(mar = c(3, 3, 3, 3))
layout(mat = matrix(c(1, 1, 2, 3, 4, 5), nrow = 2))
plot_3d_meta(meta_pm10, "pm10")
plot_3d("Newyork", "pm10")
plot_3d("Miami", "pm10")
plot_3d("Chicago", "pm10")
dev.off()

make_png("pm10", 2000, 900, "Overall")
par(mar = c(6, 5, 5, 5))
layout(mat = matrix(c(1, 1, 2, 3, 4, 5), nrow = 2))
plot_overall_meta(meta_pm10, usa_death2, "pm10", TeX("$PM_{10} \\, (\\mu g /m^3)$"))
plot_overall(usa_death2, "Newyork", "pm10", TeX("$PM_{10} \\, (\\mu g /m^3)$"))
plot_overall(usa_death2, "Miami", "pm10", TeX("$PM_{10} \\, (\\mu g /m^3)$"))
plot_overall(usa_death2, "Chicago", "pm10", TeX("$PM_{10} \\, (\\mu g /m^3)$"))
dev.off()

p1 <- plot_slice_meta(usa_death2, meta_pm10, "pm10",  "PM_{10}")
p2 <- plot_slice(usa_death2, "Newyork", "pm10",  "PM_{10}")
p3 <- plot_slice(usa_death2, "Miami", "pm10",  "PM_{10}")
p4 <- plot_slice(usa_death2, "Chicago", "pm10",  "PM_{10}")
p1 + (p2 / p3) + (p4 / plot_spacer()) + plot_layout(ncol = 3)
make_png2("pm10", 16, 8)
library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)
library(grid)
library(gridExtra)

theme_fig2 <- function() {
    theme_bw() +
        theme(
            strip.background=element_blank(),
            # strip.text=element_blank(),
            strip.text=element_text(size=8.5),
            panel.spacing.x=unit(0.1, "cm"),
            plot.margin=unit(c(0, 0, 0, 0), "in"),
            legend.title=element_text(size=7.2),
            legend.key.size=unit(4, "mm"),
            legend.text=element_text(size=7.2),
            legend.position=c(0.07, 0.83),
            legend.background=element_rect(color="gray"),
            legend.margin=margin(0.2, 0.3, 0.2, 0.3, "mm")
        )
}

texlab = c(
    "Infection rate $\\beta$"="beta", 
    "Basic rep. number $\\mathcal{R}$"="rep-number", 
    "Growth rate $\\mathcal{G}$"="growth-rate"
)

res <- read_csv("results/pi-spread-obs.csv") |>
    mutate(lab=fct_relevel(fct_recode(lab, !!!texlab), !!!names(texlab)))

res_approx <- read_csv("results/approx-pi-spread.csv") |> 
    mutate(lab=fct_relevel(fct_recode(lab, !!!texlab), !!!names(texlab)))

gg <- ggplot(res, aes(max_t, md, col=as.factor(num_obs))) +
    geom_line(data=res_approx, linetype="dashed", linewidth=1.05, alpha=0.66) +
    geom_line(alpha=0.7, linewidth=1.25)  +
    facet_wrap(~lab, scales="free") +
    scale_color_manual(values=c("lightblue", "#f53db5")) +
    coord_cartesian(ylim=c(0, NA)) +
    labs(x="Days of observation", y="Identifiability $\\delta_u$", col="Num. obs.") +
    theme_fig2()

tikz_plot(gg, "spread-obs", 4.8, 2.3)

# --------------------------------------------------------------------------------
# plot-pi-over-true-params.R------------------------------------------------------
# Produce Figure 3. Depends on some functions from "fig1-insets.R" etc------------
# --------------------------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(cowplot)
library(grid)
library(gridExtra)
library(VGAM) # for lambertW

theme_fig2 <- function() {
    theme_bw() +
        theme(
            strip.background=element_blank(),
            # strip.text=element_blank(),
            strip.text=element_text(size=8.5),
            panel.spacing.x=unit(0.1, "cm"),
            plot.margin=unit(c(0, 0, 0, 0), "in"),
            legend.position="none"
        )
}

texlab = c(
    "Recovery rate $\\alpha$"="α", 
    "Infection rate $\\beta$"="β", 
    "Initial susceptible $S_0$"="S₀", 
    "Basic rep. number $\\mathcal{R}$"="rep-number", 
    "Outbreak size $\\mathcal{O}$"="outbreak-size",
    "Peak intensity $\\mathcal{P}$"="peak-intensity",
    "Growth rate $\\mathcal{G}$"="growth-rate"
)

res <- read_csv("results/pi-over-true-by-peak.csv") |> 
    rename(alpha=α, beta=β, S0=`S₀`) |> 
    mutate(
        Reff=beta*S0/alpha,
        rnot=rep_number(alpha, beta),
        osize=outbreak_size(alpha, beta, S0),
        imax=peak_intensity(alpha, beta, S0),
        tpeak=t,
        grate=growth_rate(alpha, beta, S0),
        lab=fct_relevel(fct_recode(lab, !!!texlab), !!!names(texlab))
    ) |> 
    filter(tpeak < 29, Reff > 1)

# check which true vars were correlated with learning each var
res |> 
    group_by(lab) |> 
    summarise(across(rnot:grate, ~cor(.x, md))) |> 
    rowwise() |> 
    mutate(x=which.max(c_across(rnot:grate)))

# get correlations for Outbreak size
cor_labs <- res |> 
    group_by(lab) |> 
    summarise(cor=round(cor(osize, md), 2)) |> 
    mutate(x=-Inf, y=Inf, label=str_c("$\\rho=", cor, "$"))

gg <- res |> 
    ggplot(aes(osize, md, col=grate)) +
    # geom_line(aes(group=as.factor(S0), col=S0), alpha=0.4) +
    geom_point(alpha=0.8) +
    geom_text(
        aes(x, y, label=label), data=cor_labs, 
        inherit.aes=FALSE, size=2.8, col="gray30", hjust=-0.1, vjust=1.6
    ) +
    facet_wrap(~lab, scales="free_y", nrow=1) +
    labs(x="True outbreak size", y=NULL, col="$\\mathcal{G}^*$") +
    # scale_color_viridis_c(option="turbo") +
    scale_color_gradient(low="blue", high="#f1b3db") +
    theme_fig2() +
    theme(legend.position="none")

legend <- get_legend(
    gg + theme(
        legend.position="right",
        legend.key.height=unit(1.2, "cm"), 
        legend.key.width=unit(0.25, "cm"),
        legend.box.margin=margin(0, 0, 0, 0.3, "in")
    ))

gg <- plot_grid(gg, legend, nrow=1, rel_widths=c(1, 0.08))

y_lab <- textGrob("Identifiability $\\delta_u$", rot=90, gp=gpar(fontsize=13.5))
ggg <- grid.arrange(arrangeGrob(gg, left=y_lab))

tikz_plot(plot_grid(ggg), "pi-over-true-new", 11, 2)
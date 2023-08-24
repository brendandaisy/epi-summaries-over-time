# --------------------------------------------------------------------------------
# plot-pi-over-time.R-------------------------------------------------------------
# produce main panels of Figure 1, and combine with insets in a grid--------------
# --------------------------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)
library(grid)
library(gridExtra)

# function to produce "latexy" pdfs
tikz_plot <- function(ggp, fname='tikz', w=8.5, h=4, dir="figs") {
    cur_wd <- getwd()
    setwd(paste0(cur_wd, "/", dir))
    tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
    print(ggp)
    dev.off()
    system(paste0('lualatex ', fname, '.tex'))
    setwd(cur_wd)
}

md_plot_var <- function(var, idx) {
    res |>
        filter(var == !!var) |> 
        ggplot(aes(t, md)) +
        geom_segment(x=3, xend=3, y=-0.1, yend=1.1, col="lightblue", size=0.9, alpha=0.5, linetype="dashed") +
        geom_segment(x=8, xend=8, y=-0.1, yend=1.1, col="#f53db5", size=0.9, alpha=0.5, linetype="dashed") +
        geom_line(size=1.6, col="gray30") +
        labs(title=var, x=NULL) +
        ylim(0, 4.5) +
        theme_bw() +
        theme(
            strip.background=element_blank(),
            axis.title=element_blank(),
            axis.text=element_text(size=rel(1.1)),
            panel.grid.minor.x=element_blank(),
            legend.position="none",
            plot.margin = unit(c(t=0.03, r=0.06, b=0, l=0.09), "in")
        )
}

true_inf <- read_csv("results/inf-true.csv") |> 
    mutate(t=seq(0, 30, 0.2))

peak <- which.max(true_inf$inf)
tsteps <- seq(0, 30, 1)

texlab = c(
    "Recovery rate $\\alpha$"="α", 
    "Infection rate $\\beta$"="β", 
    "Initial susceptible $S_0$"="S₀", 
    "Basic rep. number $\\mathcal{R}$"="rep-number", 
    "Outbreak size $\\mathcal{O}$"="outbreak-size",
    "Peak intensity $\\mathcal{P}$"="peak-intensity",
    "Peak timing $\\mathcal{T}$"="peak-timing",
    "Growth rate $\\mathcal{G}$"="growth-rate"
)

res <- read_csv("results/res-pi-over-time.csv", na="missing") |> 
    mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

plot_md <- imap(names(texlab), md_plot_var)

# run a separate script to get a list of plots for the insets, then add them in corners:
source("scripts/fig1-insets.R")

md_rows <- map(1:8, ~{
    gdraw <- ggdraw(plot_md[[.x]])
    if (.x < 6)
        # gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.09, y=0.88, width=0.49, height=0.39, vjust=1)
        gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.16, y=0.88, width=0.47, height=0.37, vjust=1)
    else
        gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.96, y=0.1, width=0.47, height=0.37, hjust=1)
    # gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.9, y=0.1, width=0.49, height=0.39, hjust=1)
    gdraw
})

# arrange plots in grid and add axis labels
panel1 <- plot_grid(plotlist=md_rows, align="h", axis="l", nrow=2)
y_lab <- textGrob("Identifiability $\\delta_u$", rot=90, gp=gpar(fontsize=16))
x_lab <- textGrob("Days of observation", gp=gpar(fontsize=16))
panel1 <- grid.arrange(arrangeGrob(panel1, left=y_lab, bottom=x_lab))

tikz_plot(plot_grid(panel1), "increasing-tspan-p1", w=7.95, h=5.44)

# make part A inf. curve as a separate figure:
inf_pri <- read_csv("data/sim-full-prior.csv") |> 
    mutate(id=1:n()) |> 
    pivot_longer(-id, names_to="t", values_to="inf") |> 
    mutate(t=as.double(t))

inf_pri_summ <- inf_pri |> 
    group_by(t) |> 
    summarise(ymin=quantile(inf, 0.025), ymax=quantile(inf, 0.975))

panel2 <- ggplot(true_inf, aes(t, inf)) +
    annotate("rect", xmin=-1, ymin=0, xmax=6, ymax=0.5, fill="lightblue", alpha=0.35) +
    annotate("rect", xmin=6, ymin=0, xmax=14, ymax=0.5, fill="#f53db5", alpha=0.2) +
    annotate("rect", xmin=14, ymin=0, xmax=31, ymax=0.5, fill="blue4", alpha=0.2) +
    geom_line(col="orange", size=1.4) +
    geom_point(data=tibble(t=0:30, inf=stan_dat$y/1000), col="#f53db5", size=0.9) +
    scale_x_continuous(guide=guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    labs(x="Days $t$", y="$I(t)$", col=NULL) +
    theme_half_open() +
    theme(plot.margin=unit(c(0, 0, 0, 0.25), "in"))

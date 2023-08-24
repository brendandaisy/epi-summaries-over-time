# --------------------------------------------------------------------------------
# plot-pi-spread-obs.R------------------------------------------------------------
# short script to make Figure 2, evenly spreading obs. over diff. timespans-------
# --------------------------------------------------------------------------------

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
            plot.margin=unit(c(0, 0.05, 0, 0.1), "in"),
            legend.title=element_text(size=7.2),
            legend.key.size=unit(4, "mm"),
            legend.text=element_text(size=7.2),
            legend.position=c(0.21, 0.85),
            legend.background=element_rect(color="gray"),
            legend.margin=margin(0.2, 0.3, 0.2, 0.3, "mm"),
            plot.title=element_text(size=9.5)
        )
}

spread_obs_facet <- function(var, pri_dens, breaks=waiver()) {
    res_sub <- filter(res, lab == var)
    res_apx_sub <- filter(res_approx, lab == var)
    
    ggplot(res_sub, aes(max_t, md, col=as.factor(num_obs))) +
        geom_line(data=res_apx_sub, linetype="dashed", linewidth=1.05, alpha=0.66) +
        geom_line(alpha=0.7, linewidth=1.25)  +
        scale_color_manual(values=c("lightblue", "#f53db5")) +
        scale_y_continuous(sec.axis=sec_axis(~ 1 / (exp(.) * sqrt(2*pi) * pri_dens), breaks=breaks)) +
        coord_cartesian(ylim=c(0, NA)) +
        labs(x=NULL, y=NULL, col="Num. obs.", title=var) +
        theme_fig2()
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

ppp <- plot_grid(
    spread_obs_facet("Infection rate $\\beta$", dunif(1.25, 0.3, 1.5)),
    spread_obs_facet("Basic rep. number $\\mathcal{R}$", 0.03, round(seq(14, 1, length.out=5), 1)) + theme(legend.position="none"),
    spread_obs_facet("Growth rate $\\mathcal{G}$", 0.41, c(round(seq(1.1, 0.01, length.out=5)[1:4], 1), 0.02)) + theme(legend.position="none"),
    nrow=1
)

y_lab1 <- textGrob("Identifiability $\\delta_u$", rot=90, gp=gpar(fontsize=11))
y_lab2 <- textGrob("Approx. std. err.", rot=270, gp=gpar(fontsize=11))
x_lab <- textGrob("Span of observation (days)", gp=gpar(fontsize=11))
pppp <- grid.arrange(arrangeGrob(ppp, left=y_lab1, bottom=x_lab, right=y_lab2))

tikz_plot(plot_grid(pppp), "spread-obs-new", 6.2, 2.45)

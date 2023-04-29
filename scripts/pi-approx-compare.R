library(tidyverse)
library(ggthemes)
library(cowplot)
library(tikzDevice)
library(grid)
library(gridExtra)

tikz_plot <- function(ggp, fname='tikz', w=8.5, h=4, dir="figs") {
    cur_wd <- getwd()
    setwd(paste0(cur_wd, "/", dir))
    tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
    print(ggp)
    dev.off()
    system(paste0('lualatex ', fname, '.tex'))
    setwd(cur_wd)
}

res <- read_csv("results/pi-approx-ts.csv") |> 
    # mutate(nrep=c(1, 5, 10)) |> 
    pivot_longer(-ts)

res_approx <- read_csv("results/approx-pi-curve-ts.csv") |> 
    mutate(ts=seq(0.2, 2, 0.1))

res_S0 <- filter(res, name == "mdS0")
res_grate <- filter(res, name == "mdgrate")

gg1 <- ggplot(res_S0, aes(1/ts, value)) +
    geom_line(aes(y=mdS0_approx), res_approx, col="#a095d5", linewidth=1.3) +
    geom_point(size=1.2, col="gray30", shape=2) +
    # stat_function(fun=~0.5 * log(.x) + res$value[3], linetype="dashed", col="gray70", linewidth=1.1) +
    labs(y="$\\delta_{S_0}$", x=NULL) +
    theme_bw()

gg2 <- ggplot(res_grate, aes(1/ts, value)) +
    geom_line(aes(y=mdgrate_approx), res_approx, col="#a095d5", linewidth=1.3) +
    geom_point(size=1.2, col="gray30", shape=2) +
    # stat_function(fun=~0.5 * log(.x) + res$value[4], linetype="dashed", col="gray70", linewidth=1.1) +
    labs(y="$\\delta_{\\mathcal{G}}$", x=NULL) +
    theme_bw()

gg <- plot_grid(gg1, gg2)

x_lab <- textGrob("Observation rate (days$^{-1}$)", gp=gpar(fontsize=11))
ggg <- grid.arrange(arrangeGrob(gg, bottom=x_lab))

tikz_plot(plot_grid(ggg), "pi-approx-compare", 4.2, 2.4)

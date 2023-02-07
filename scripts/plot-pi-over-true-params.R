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

res <- read_csv("data/pi-over-true-by-peak.csv") |> 
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

# this of course reflects association controlling for effect of other variables
# TODO also these are very correlated so these estimates are meaningless anyways!!!!!!!!!
res_lm <- res |> 
    group_by(lab) |> 
    group_modify(~{
        fit <- lm(md ~ 0 + rnot + osize + imax + tpeak + grate, data=.x)
        as_tibble(summary(fit)$coefficients, rownames="predictor")
    }, .keep=TRUE) |> 
    ungroup()

res_lm |> 
    group_by(predictor) |> 
    summarize(num_sig=sum(`Pr(>|t|)` < 0.05), num_pos=sum(Estimate < 0))

res |> 
    group_by(lab) |> 
    summarise(across(rnot:grate, ~cor(.x, md))) |> 
    rowwise() |> 
    mutate(x=which.max(c_across(rnot:grate)))

cor_labs1 <- res |> 
    group_by(lab) |> 
    summarise(cor=round(cor(tpeak, md), 2)) |> 
    mutate(x=Inf, y=Inf, label=str_c("$\\rho=", cor, "$"))

cor_labs <- res |> 
    group_by(lab) |> 
    summarise(cor=round(cor(osize, md), 2)) |> 
    mutate(x=-Inf, y=Inf, label=str_c("$\\rho=", cor, "$"))

ggplot(res, aes(osize, tpeak))

p1 <- ggplot(res, aes(tpeak, md, col=beta)) +
    # geom_line(aes(group=as.factor(S0), col=S0), alpha=0.4) +
    geom_point() +
    geom_text(
        aes(x, y, label=label), data=cor_labs1, 
        inherit.aes=FALSE, size=2.8, col="#dc4a29", hjust=1.1, vjust=1.6
    ) +
    facet_wrap(~lab, scales="free_y", nrow=1) +
    labs(x="Peak timing $\\mathcal{T}^*$", y=NULL, col="$\\beta^*$") +
    theme_fig2()

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
    theme(legend.position="right")

legend <- get_legend(
    p1 + theme(
        legend.position="right",
        legend.key.height=unit(1.2, "cm"), 
        legend.key.width=unit(0.25, "cm"),
        legend.box.margin=margin(0, 0, 0, 0.3, "in")
    ))

gg <- plot_grid(
    plot_grid(p1, p2, align="h", nrow=2),
    legend,
    nrow=1,
    rel_widths=c(1, 0.08)
)

y_lab <- textGrob("Identifiability $\\delta_u$", rot=90, gp=gpar(fontsize=13.5))
ggg <- grid.arrange(arrangeGrob(gg, left=y_lab))

tikz_plot(plot_grid(ggg), "pi-over-true-new", 11, 2)

###

ml_test <- read_csv("data/ml-bench-res.csv")
    
p1 <- ggplot(ml_test, aes(M, ((10^-9)*mean_time))) +
    geom_line(col="#a095d5", linewidth=1.4) +
    # facet_wrap(~name, scales="free") +
    theme_bw() +
    scale_x_continuous(labels=scales::label_scientific(digits=1)) +
    labs(y="Avg. runtime (sec)", x=NULL) +
    theme(plot.margin=unit(c(0.2, .5, 0, 0), "cm"))

p2 <- ggplot(ml_test, aes(M, std_res)) +
    geom_line(col="#a095d5", linewidth=1.4) +
    # facet_wrap(~name, scales="free") +
    theme_bw() +
    scale_x_continuous(labels=scales::label_scientific(digits=1)) +
    # coord_trans(y="log10") +
    scale_y_sqrt(breaks=c(1, seq(10, 60, 10))) +
    # annotation_logticks(scaled=TRUE, sides="l") +
    labs(y="Std. error", x=NULL) +
    theme(plot.margin=unit(c(0.2, 0.3, 0, 0.1), "cm"))

gg <- plot_grid(p1, p2, nrow=1)
x_lab <- textGrob("Monte Carlo iterations", gp=gpar(fontsize=10))
ggg <- grid.arrange(arrangeGrob(gg, bottom=x_lab))

tikz_plot(plot_grid(ggg), "marg-lik-std-err", w=4, h=2.2)

###

res_sub <- filter(res, str_detect(lab, "alpha"))
lfit <- lm(md ~ tpeak + osize, data=res_sub)

# predict over sensible grid of values
tpeak <- unique(res$tpeak)
osize <- unique(res$osize)

pred_grid <- expand_grid(tpeak, osize)
pred_mat <- matrix(predict(lfit, newdata=pred_grid), nrow=length(tpeak), ncol=length(osize))


plot_ly() |> 
    add_surface(x=tpeak, y=osize, z=pred_mat) |> 
    add_trace(x=~osize, y=~tpeak, z=~md, color=~beta, data=res_sub, type="scatter3d")

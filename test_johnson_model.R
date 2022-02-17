# Test Johnson model

library(ggplot2)
library(ggthemes)
library(patchwork)

source("johnson_model.R")

# run static model
par_static <- model(PAR = seq(1, 2400, 24), CB6F = (1.2/1e6), RUB = (27.8/1e6))

# run dynamic model
par_dynamic <- model(PAR = seq(1, 2400, 24), CB6F = (1.2/1e6), RUB = (27.8/1e6), Ku2 = 2e09, 
                     alpha_opt = "dynamic")

# plot static model outputs
(Photosystems_static <- ggplot() + 
    geom_line(data = par_static, aes(x = PAR, y = a2/Abs*100, color = "Photosystem II")) +
    geom_line(data = par_static, aes(x = PAR, y = a1/Abs*100, color = "Photosystem I")) +
    theme_few() +
    scale_color_manual(name = NULL, breaks = c("Photosystem II", "Photosystem I"),
                       values = c("Photosystem II" = "blue", "Photosystem I" = "red")) +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Cross-sections (%)", limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)))
(PAM1_static <- ggplot() + 
    geom_line(data = par_static, aes(x = PAR, y = PAM1_a)) +
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Photochemical Yield", limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
(PAM2_static <- ggplot() + 
    geom_line(data = par_static, aes(x = PAR, y = PAM2_a)) +
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Reg. NPQ Yield", limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
(PAM3_static <- ggplot() + 
    geom_line(data = par_static, aes(x = PAR, y = PAM3_a)) +    
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Residual Yield", limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
(A_static <- ggplot() + 
    geom_line(data = par_static, aes(x = PAR, y = An_a * 1e6)) +
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "An (umol CO2 m-2 s-1", limits = c(0, 30), breaks = c(0, 10, 20, 30)))

(static_figures <- Photosystems_static + PAM1_static + PAM2_static + PAM3_static + A_static + 
    plot_layout(ncol = 1, nrow = 5))

# plot dynamic model outputs
(Photosystems_dynamic <- ggplot() + 
    geom_line(data = par_dynamic, aes(x = PAR, y = a2/Abs*100, color = "Photosystem II")) +
    geom_line(data = par_dynamic, aes(x = PAR, y = a1/Abs*100, color = "Photosystem I")) +
    theme_few() +
    scale_color_manual(name = NULL, breaks = c("Photosystem II", "Photosystem I"),
                       values = c("Photosystem II" = "blue", "Photosystem I" = "red")) +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Cross-sections (%)", limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)))
(PAM1_dynamic <- ggplot() + 
    geom_line(data = par_dynamic, aes(x = PAR, y = PAM1_a)) +
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Photochemical Yield", limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
(PAM2_dynamic <- ggplot() + 
    geom_line(data = par_dynamic, aes(x = PAR, y = PAM2_a)) +
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Reg. NPQ Yield", limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
(PAM3_dynamic <- ggplot() + 
    geom_line(data = par_dynamic, aes(x = PAR, y = PAM3_a)) +    
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "Residual Yield", limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
(A_dynamic <- ggplot() + 
    geom_line(data = par_dynamic, aes(x = PAR, y = An_a * 1e6)) +
    theme_few() +
    scale_x_continuous(name = "PAR (umol PPFD m-2 s-1)", breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
    scale_y_continuous(name = "An (umol CO2 m-2 s-1", limits = c(0, 30), breaks = c(0, 10, 20, 30)))

(dynamic_figures <- Photosystems_dynamic + PAM1_dynamic + PAM2_dynamic + PAM3_dynamic + A_dynamic + 
    plot_layout(ncol = 1, nrow = 5))

# save figures
ggsave("example_figures/figures_static.jpeg", static_figures, 
       width = 500, height = 1000, units = "px", dpi = 72)
ggsave("example_figures/figures_dynamic.jpeg", dynamic_figures, 
       width = 500, height = 1000, units = "px", dpi = 72)

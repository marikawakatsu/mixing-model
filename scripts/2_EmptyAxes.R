################################################################################
#
# Make empty plot for Fig. 4
#
################################################################################

rm(list = ls())

source("scripts/util/__Util__MASTER.R")

shape <- 4
colour = "black"

# emplot <- ggplot() +
#   xlab( expression("Demand rate"~(delta) ) ) +
#   ylab( expression( Delta~"mean threshold"~(mu^X - mu^Y) ) ) +
#   theme_mk() +
#   # theme(axis.line = element_line(colour = 'black' )) +
#   scale_x_continuous(limits = c(0.3, 1.8), breaks = c(0.6, 1.5)) +
#   scale_y_continuous(limits = c(-2, 12), breaks = c(0, 4, 10)) +
#   geom_point(aes(x = 0.6, y = 0), colour = colour, shape = shape) +
#   geom_point(aes(x = 0.6, y = 4), colour = colour, shape = shape) +
#   geom_point(aes(x = 0.6, y = 10), colour = colour, shape = shape) +
#   geom_point(aes(x = 1.5, y = 0), colour = colour, shape = shape)

emplot <- ggplot() +
  xlab( expression("Demand rate"~(delta) ) ) +
  ylab( expression( Delta~"mean threshold"~(mu^Y - mu^X) ) ) +
  theme_mk() +
  # theme(axis.line = element_line(colour = 'black' )) +
  scale_x_continuous(limits = c(0.3, 1.8), breaks = c(0.6, 1.5)) +
  scale_y_continuous(limits = c(-12, 2), breaks = c(-10, -4, 0)) +
  geom_point(aes(x = 0.6, y = 0), colour = colour, shape = shape) +
  geom_point(aes(x = 0.6, y = -4), colour = colour, shape = shape) +
  geom_point(aes(x = 0.6, y = -10), colour = colour, shape = shape) +
  geom_point(aes(x = 1.5, y = 0), colour = colour, shape = shape)

emplot

ggsave(filename = paste0("slides/Fig4_EmptyAxes.pdf"), width = 1.6, height = 1.5, dpi = 800)




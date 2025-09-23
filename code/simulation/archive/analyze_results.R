library(here)
library(data.table)
library(ggplot2)

source(here("R/population_obj.R"))
source(here("R/sim_data.R"))

df <- lapply(list.files("outputs/01_highX1_linear/", full.names = TRUE), function(file){
  dt <- readRDS(file)
  return(dt$est)
})

df <- lapply(list.files("outputs/02_equal_linear/", full.names = TRUE), function(file){
  dt <- readRDS(file)
  return(dt$est)
})

dt <- rbindlist(df)


plot_dt <- dt[method %in% c("cs", "pCCA")]

#--- new vectors to add -------------------------------------------------
# ensure plain numeric length-2 vectors
sim        <- sim_linear(n = 200,
                         Sigma_e = diag(c(9, 1)),
                         beta = c(1,0.1)
)
omega <- sim$omega
omega.do   <- optim_dopca(beta, sim$Sigma_e,sim$Sigma_X)
omega.do.s <- optim_dopca(beta, sim$Sigma_e,sim$Sigma_X, TRUE)

add_dt <- data.table(
  method = c("omega.cs.true", "omega.do.opt", "omega.do.s.opt"),
  X1     = c(omega.cs[1],   omega.do[1],   omega.do.s[1]),
  X2     = c(omega.cs[2],   omega.do[2],   omega.do.s[2])
)

# optional: factor levels so legend orders nicely
plot_dt[, method := factor(method, levels = c("cs","do","do_scaled","pCCA"))]
add_dt[,  method := factor(method, 
                           levels = c(levels(plot_dt$method), 
                                      "omega.cs.true","omega.do.opt","omega.do.s.opt"))]

#--- build plot ---------------------------------------------------------
p <- ggplot(plot_dt, aes(x = 0, y = 0, xend = X1, yend = X2, col = method)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), alpha = 0.3) +
  # add the 3 highlighted vectors
  geom_segment(data = add_dt,
               aes(x = 0, y = 0, xend = X1, yend = X2, col = method),
               arrow = arrow(length = unit(0.25, "cm")),
               size = 1.1, alpha = 1) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Estimated Vectors and True / DO Optima",
       x = "X", y = "Y") +
  theme(legend.title = element_blank())

p

ggsave("estimated_omegas.png")

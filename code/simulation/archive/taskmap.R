library(data.table)

B <- 100

taskmap <- CJ(
  n = c(200, 800, 1600),
  scenario = 1:5,
  i = 1:B
)

# Add interaction model at the end....keep it linear though. 
# So model 1 with interaction. 

taskmap2 <- 

taskmap[, ID := 1:.N]




fwrite(taskmap, "config/taskmap.csv")

dyn.load('./C/bekk_log_lik.so')
source('./R/bekk_model.R')

y <- read.csv('../ganter/data/Problemsets Data (Tickers)/all.csv')
y <- apply(y,c(2),diff)
y[is.na(y)] <- 0

fit <- scalar.bekk.fit(as.matrix(y[,2:4]),opts=list(fit = TRUE))
print(fit$param)
print(fit$obj)

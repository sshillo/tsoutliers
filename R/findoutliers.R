# library(forecast)
# library(tsoutliers)
# library(xts)
# df <- read.csv("/Users/sshillo/Development/statsketch/backend/statsketch/fixtures/TLC_rej_daily_summary.csv")
# df$date = as.Date(df$date,format="%m/%d/%y")
# xdf <- xts(df$walkby_conv,order.by=df$date)
# tdf <- as.ts(xdf)
# fit <- auto.arima(tdf)
# resid <- residuals(fit)
# pars <- coefs2poly(fit)
# outliers <- locate.outliers(resid, pars,10)
# filtered <- outliers[order(-abs(outliers$tstat)),][1:5,]
# 

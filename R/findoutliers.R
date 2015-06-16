#library(forecast)
#library(tsoutliers)
#library(xts)
#df <- read.csv("/Users/sshillo/Development/statsketch/backend/statsketch/fixtures/TLC_rej_daily_summary.csv")
#df$app_date = as.Date(df$app_date)
#xdf <- xts(df$avg_debt_to_income,order.by=df$app_date)
#tdf <- as.ts(xdf)
#fit <- auto.arima(tdf)
#resid <- residuals(fit)
#pars <- coefs2poly(fit)
#outliers <- locate.outliers(resid, pars,10)
#filtered <- outliers[order(-abs(outliers$tstat)),][1:5,]


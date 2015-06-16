TsoFastBase <- function(x, xreg = NULL, cval = 3.5, delta = 0.7, n.start = 50,
                 types = c("AO", "LS", "TC"), maxit.iloop = 4, 
                 remove.method = c("en-masse", "bottom-up"),
                 remove.cval = NULL, 
                 tsmethod = c("auto.arima", "arima", "stsm"), args.tsmethod = NULL,
                 args.tsmodel = NULL, logfile = NULL)
{
  # "x" can be either a "ts" object or a "stsm" object;
  # if !inherits(x, "stsm") then two identical objects are stored ("x" and "y")
  y <- if(is.ts(x)) { x } else x@y
  
  #remove.method <- match.arg(remove.method)
  #tsmethod <- match.arg(tsmethod)
  #remove.method <- match.arg(remove.method)
  fitmethod <- gsub("stsm", "stsmFit", tsmethod)
  
  if (is.null(remove.cval))
    remove.cval <- cval
  
  # fit time series model
  
  fit <- do.call(fitmethod, args = c(list(x = x, xreg = xreg), args.tsmethod))
  #fit$series <- deparse(substitute(y))
  
  
  if (!is.null(logfile))
  {
    cat(paste("model selection:\n"), file = logfile, append = FALSE)
    capture.output(fit, file = logfile, append = TRUE)
  }
  
  # identify and locate prospective outliers by type
  # given a fitted time series model
  stage1 <- locate.outliers.oloop(y = y, fit = fit, types = types, cval = cval, 
                                  maxit.iloop = maxit.iloop, delta = delta, n.start = n.start, logfile = logfile)
  
  
  # choose and fit the model including the outlier regressors detected so far
  # (the weights of the outliers is fine tuned, to see it 
  # compare 'moall[,"coefhat"]' with 'coef(fit)["oeffi"]') then
  # remove the outliers detected so far if they are not significant in the new model/fit
  
  if (nrow(stage1$outliers) > 0)
  {
    stage1$outliers = stage1$outliers[order(-abs(stage1$outliers$tstat)),][1:20,]
    stage2 <- remove.outliers(x = stage1, y = y, cval = remove.cval, 
                              method = remove.method, delta = delta, n.start = n.start, 
                              tsmethod.call = fit$call, fdiff = NULL, logfile = logfile)
    
    stopifnot(ncol(stage2$xreg) == length(stage2$xregcoefs))
  } else 
    stage2 <- list(xreg = NULL, fit = stage1$fit)
  
  # final outliers and
  # linearized series, original series without outlier effects
  
  if (!is.null(stage2$xreg))
  {
    # stage2$fit$xreg is not returned by arima()
    moall <- stage2$outliers
    moall[,"coefhat"] <- stage2$xregcoefs
    moall[,"tstat"] <- stage2$xregtstats
    
    oeff <- stage2$xreg %*% cbind(stage2$xregcoefs)
    attributes(oeff) <- attributes(y)
    yadj <- y - oeff
    
    moall <- moall[,c("type", "ind", "coefhat", "tstat")]
    outtimes <- time(y)[moall[,"ind"]]
    if (frequency(y) > 1) 
      outseason <- formatC(as.vector(cycle(y)[moall[,"ind"]]), 
                           width = 2, flag="0")
    
    moall <- cbind(moall[,c("type", "ind")], 
                   "time" = if (frequency(y) > 1) paste(floor(outtimes), 
                                                        outseason, sep = ":") else outtimes,
                   moall[,c("coefhat","tstat")])
    
    oind <- order(moall[,"ind"])
    moall <- moall[oind,]
    outtimes <- outtimes[oind]
    rownames(moall) <- NULL
    
  } else { # no outliers detected
    oeff <- NULL
    yadj <- y
    moall <- data.frame(array(dim = c(0, 4)))
    colnames(moall) <- c("type", "ind", "coefhat", "tstat")
    outtimes <- NULL
  }
  
  if (!is.null(logfile))
  {
    msg <- paste("\nfinal outliers\n")
    cat(msg, file = logfile, append = TRUE)
    capture.output(moall, file = logfile, append = TRUE)
  }
  structure(list(outliers = moall, y = y, yadj = yadj, cval = cval,
                 fit = stage2$fit, effects = oeff, times = outtimes), 
            class = "tsoutliers")
}

plot.cumres <- function(x, idx=1:length(x$variable),
                        col=c("grey"),
                        ci=TRUE, col.ci="darkblue", col.alpha=0.3, lty.ci=0, level=0.95,
                        legend=c("type1","type2","none"),
                        ylim=NULL, ...) {
  ylim. <- ylim
  for (i in idx) {
    legendtxt <- c(); legendpch <- c(); legendcol <- c(); legendlty <- c(); legendlwd <- c(); legendcex <- c()
    if (is.null(ylim)) {
      ylim. <- max(abs(range(x$W[,i])))*2*c(-1,1)
    }
    ## Observed process
    xlab <- x$variable[i]; main <- ""
    if (x$type[i]=="score") {
      main <- xlab; xlab <- "Time";
    }
    legendtxt <- c(legendtxt, "Observed"); legendpch <- c(legendpch,-1); legendcol <- c(legendcol,1); legendlty <- c(legendlty,1); legendlwd <- c(legendlwd,2); legendcex <- c(legendcex,1);
    with(x, plot(W[,i] ~ x[,i], type="s", lwd=2, ylab=substitute(expression(W[p](x)),list(p=variable[i])),ylim=ylim.,xlab=xlab,main=main));

    ## Sample processes
    if (col!="none" && !is.null(col)) {
      legendtxt <- c(legendtxt, "MC sample"); legendpch <- c(legendpch,-1); legendcol <- c(legendcol,col); legendlty <- c(legendlty,1); legendlwd <- c(legendlwd,1); legendcex <- c(legendcex,1);
      for (k in 1:ncol(x$What[[i]])) {
        lines(x$What[[i]][,k] ~ x$x[,i], type="s", col=col, lwd=1)
      }; lines(x$W[,i] ~ x$x[,i], type="s", lwd=2)
    }

    ## Prediction bandds
    if ( (ci[1]!="none" && !is.null(ci) && ci[1]!=0) || (ci==TRUE) ) {
##    if ((ci[1]!="none" && !is.null(ci))) {
      if (col.alpha==0)
        col.trans <- col.ci
      else 
        col.trans <- sapply(col.ci, FUN=function(x) do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
      if (ci[1]=="pointwise")
        myCI <- confint(x,parm=i,cval=qnorm(1-(1-level)/2))
      else
        myCI <- confint(x,parm=i,level=level)      
      mystepf <- with(myCI, stepfun(t,c(0,yu)));
      t <- c();
      epsilon <- 1e-9
      for (k in 1:length(myCI$t)) {
        t <- c(t, myCI$t[k]-epsilon, myCI$t[k])
      }; t <- t[-1]
      yu <- mystepf(t)
      legendtxt <- c(legendtxt, "95% prediction band"); legendpch <- c(legendpch,15); legendcol <- c(legendcol,col.trans); legendlty <- c(legendlty,0); legendlwd <- c(legendlwd,0); legendcex <- c(legendcex,2);
      
      lines(yu ~ t, lwd=1, col=col.ci, lty=lty.ci)
      lines(-yu ~ t, lwd=1, col=col.ci, lty=lty.ci)
      tt <- c(t, rev(t))
      yy <- c(yu, rev(-yu))
      polygon(tt,yy, col=col.trans, lty=0)      
      ##as.list(c(col2rgb("darkblue"),10)/255))),
    }
    
    if (!is.null(legend) && legend[1]!="none") {
      if (legend[1]=="type1")
        legend("topright", c(paste("KS-test: p=",x$KS[i],sep=""),paste("CvM-test: p=",x$CvM[i],sep="")), bg="white")
      else
        legend("topright", legendtxt, lty=legendlty, pch=legendpch, col=legendcol, lwd=legendlwd, pt.cex=legendcex, bg="white")
    }
    ylim. <- NULL
  }
  invisible()
}




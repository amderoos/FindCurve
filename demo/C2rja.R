devAskNewPage(ask = FALSE)
oldwd = getwd()
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)
graphics.off()

setwd(system.file("Systems", package="FindCurve"))

init = c(2.07278E-01, 2.07266E-01, 2.00000E+00, 1.00656E-05, 8.00686E-08);
parameters = c(1.0, 2.0, 2.0, 0.5, 0.1, 10.0, 1.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 10.0);

cat('\n\n\nConsumer-resource equilibrium as a function of R1max for default parameters\n\n');
cmd = 'output1 <- FindCurve("C2rja", "EQ", init, 0.1, c(0.0, 3.0), parameters, options = c("par1", "1"), clean = TRUE, force = TRUE)'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))
  cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

  par(mar = c(4, 5, 2, 5), tcl=0.4)
  plot(1, 1, type="l", xaxs="i", yaxs="i", xlim=c(0.0, 3.0), ylim=c(0.0, 1.2), xlab="", ylab="")
  mtext("Maximum density resource 1", 1, line=2.5, cex=1.3)
  mtext("Average consumer density", 2, line=3.5, cex=1.3)

  lines(output1$curvepoints[,1], output1$curvepoints[,9], type="l", col=rgb(0.6,0,0), lwd=3)
  points(output1$bifpoints[,1], output1$bifpoints[,9], col="red", pch=8, lwd=2)
  text(output1$bifpoints[,1], output1$bifpoints[,9], output1$biftype, pos=3, offset=0.35, adj = c(0.1))

  lines(output1$curvepoints[,1], output1$curvepoints[,10], type="l", col=rgb(0,0,0.6), lwd=3)
  points(output1$bifpoints[,1], output1$bifpoints[,10], col="red", pch=8, lwd=2)
  text(output1$bifpoints[,1], output1$bifpoints[,10], output1$biftype, pos=3, offset=0.35, adj = c(0.1))
}

par(par.defaults)
CleanFindCurve('F')
setwd(oldwd)

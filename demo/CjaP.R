devAskNewPage(ask = FALSE)
oldwd = getwd()
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)
graphics.off()

setwd(system.file("Systems", package="FindCurve"))

parameters = c(0.1, 30.0, 1.0, 0.1, 0.5, 0.1, 0.015, 0.5, 3.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.32, 0.032, 0.5, 0.005, 0.0, 3.0, 0.0, 1.0, 0.0);

cat('\n\n\nTesting juvenile biomass overcompensation detection in the consumer-resource equilibrium\n\n');
cmd = 'output1 <- FindCurve("CjaP", "EQ", c(0.0, 2.12019, 0.892170, 11.3152, 0.0), 0.5, c(0.0, 1.0), parameters, c("par1", "9", "par2", "13", "EXTfun", "1", "EXTpar", "9"), clean = TRUE, force = TRUE)'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))
  cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

  par(mar = c(4, 5, 2, 5), tcl=0.4)
  plot(1, 1, type="l", xaxs="i", yaxs="i", xlim=c(0.0, 0.3), ylim=c(0.0, 3.0), xlab="", ylab="")
  mtext("Additional consumer mortality", 1, line=2.5, cex=1.3)
  mtext("Average juvenile density", 2, line=3.5, cex=1.3)

  lines(output1$curvepoints[,1], output1$curvepoints[,8], type="l", col=rgb(0,.6,0), lwd=3)
  points(output1$bifpoints[,1], output1$bifpoints[,8], col="red", pch=8, lwd=2)
  text(output1$bifpoints[,1], output1$bifpoints[,8], output1$biftype, pos=3, offset=0.35, adj = c(0.1))
} else {
    par(par.defaults)
    CleanFindCurve('F')
    graphics.off()
    setwd(oldwd)
    stop('Test script interrupted by user')
  }

cat('\n\n\nContinuation of the juvenile biomass overcompensation boundary in the consumer-resource equilibrium as a function of consumer mortality and interval\n\n');
cmd = 'output2a <- FindCurve("CjaP", "EXT", output1$bifpoints[output1$biftypes == "EXT",1:6], 0.5, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "9", "par2", "13", "EXTfun", "1", "EXTpar", "9"))'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))

  cmd = 'output2b <- FindCurve("CjaP", "EXT", output1$bifpoints[output1$biftypes == "EXT",1:6], -0.5, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "9", "par2", "13", "EXTfun", "1", "EXTpar", "9"))'
  eval(parse(text=cmd))

  plot(1, 1, type="l", xaxs="i", yaxs="i", xlim=c(0.0, 0.3), ylim=c(0.0, 30.0), xlab="", ylab="")
  mtext("Additional consumer mortality", 1, line=2.5, cex=1.3)
  mtext("Interval", 2, line=3.5, cex=1.3)

  lines(output2a$curvepoints[,1], output2a$curvepoints[,6], type="l", col=rgb(0,.6,0), lwd=3)
  lines(output2b$curvepoints[,1], output2b$curvepoints[,6], type="l", col=rgb(0,.6,0), lwd=3)
} else {
    par(par.defaults)
    CleanFindCurve('F')
    graphics.off()
    setwd(oldwd)
    stop('Test script interrupted by user')
  }

cat('\n\n\nContinuation of the the predator extinction boundary as a function of consumer mortality and interval\n\n');
cmd = 'output3a <- FindCurve("CjaP", "BP", output1$bifpoints[output1$biftypes == "BP",1:6], 0.5, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "9", "par2", "13"))'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))

  cmd = 'output3b <- FindCurve("CjaP", "BP", output1$bifpoints[output1$biftypes == "BP",1:6], -0.5, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "9", "par2", "13"))'
  eval(parse(text=cmd))

  lines(output3a$curvepoints[,1], output3a$curvepoints[,6], type="l", col=rgb(.6,0,0), lwd=3)
  lines(output3b$curvepoints[,1], output3b$curvepoints[,6], type="l", col=rgb(.6,0,0), lwd=3)
} else {
    par(par.defaults)
    CleanFindCurve('F')
    graphics.off()
    setwd(oldwd)
    stop('Test script interrupted by user')
  }

parameters = c(0.1, 30.0, 1.0, 0.1, 0.5, 0.1, 0.015, 1.5, 3.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.32, 0.032, 0.5, 0.005, 0.0, 3.0, 0.0, 1.0, 0.0)

cat('\n\n\nTesting saddle-node bifurcation detection in the predator consumer-resource equilibrium\n\n');
cmd = 'output1 <- FindCurve("CjaP", "EQ", c(0.0, 27.0396, 0.590107, 0.132521, 1.84090), 0.2, c(0.0, 1.0), parameters, c("par1", "21", "par2", "13"))'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))
  cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

  par(mar = c(4, 5, 2, 5), tcl=0.4)
  plot(1, 1, type="l", xaxs="i", yaxs="i", xlim=c(0.0, 0.1), ylim=c(0.0, 5.0), xlab="", ylab="")
  mtext("Additional predator mortality", 1, line=2.5, cex=1.3)
  mtext("Average predator density", 2, line=3.5, cex=1.3)

  lines(output1$curvepoints[,1], output1$curvepoints[,11], type="l", col=rgb(0,.6,0), lwd=3)
  points(output1$bifpoints[,1], output1$bifpoints[,11], col="red", pch=8, lwd=2)
  text(output1$bifpoints[,1], output1$bifpoints[,11], output1$biftype, pos=3, offset=0.35, adj = c(0.1))
} else {
    par(par.defaults)
    CleanFindCurve('F')
    graphics.off()
    setwd(oldwd)
    stop('Test script interrupted by user')
  }

cat('\n\n\nContinuation of the saddle-node bifurcation as a function of predator mortality and interval\n\n');
cmd = 'output2 <- FindCurve("CjaP", "LP", output1$bifpoints[1,1:6], -0.2, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "21", "par2", "13"))'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))

  cmd = 'output3 <- FindCurve("CjaP", "LP", output1$bifpoints[2,1:6], -0.2, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "21", "par2", "13"))'
  eval(parse(text=cmd))

  plot(1, 1, type="l", xaxs="i", yaxs="i", xlim=c(0.0, 0.1), ylim=c(0.0, 70.0), xlab="", ylab="")
  mtext("Additional predator mortality", 1, line=2.5, cex=1.3)
  mtext("Interval", 2, line=3.5, cex=1.3)

  lines(output2$curvepoints[,1], output2$curvepoints[,6], type="l", col=rgb(.6,0,0), lwd=3)
  lines(output3$curvepoints[,1], output3$curvepoints[,6], type="l", col=rgb(.6,0,0), lwd=3)
} else {
    par(par.defaults)
    CleanFindCurve('F')
    graphics.off()
    setwd(oldwd)
    stop('Test script interrupted by user')
  }

cat('\n\n\nContinuation of the predator invasion boundary as a function of predator mortality and interval\n\n');
cmd = 'output4 <- FindCurve("CjaP", "BP", output1$bifpoints[3,1:6], -0.2, c(0.0, 1.0, 0.0, 70.0), parameters, c("par1", "21", "par2", "13"))'

str = readline(paste0("\n> ", cmd, " [y(es)/q(uit)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))

  lines(output4$curvepoints[,1], output4$curvepoints[,6], type="l", col=rgb(0,0,0.6), lwd=3)
} else {
    par(par.defaults)
    CleanFindCurve('F')
    graphics.off()
    setwd(oldwd)
    stop('Test script interrupted by user')
  }

par(par.defaults)
CleanFindCurve('F')
graphics.off()
setwd(oldwd)

oldwd = getwd()

parameters = c(0.1, 30.0, 1.0, 0.1, 0.5, 0.1, 0.015, 0.5, 3.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.32, 0.032, 0.5, 0.005, 0.0, 3.0, 0.0, 1.0, 0.0);

cat('\n\n\nTesting juvenile biomass overcompensation detection in the consumer-resource equilibrium\n\n');
cmd = 'output1 <- FindCurve("CjaP", "EQ", c(1.0E-4, 2.12019, 0.892170, 11.3152, 0.0), 0.5, c(1.0E-4, 1.0), parameters, c("par1", "9", "EXTfun", "1", "EXTpar", "9", "silent"), clean = TRUE, force = TRUE)'

eval(parse(text=cmd))
cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

cat('\n\n\nContinuation of the juvenile biomass overcompensation boundary in the consumer-resource equilibrium as a function of consumer mortality and interval\n\n');
cmd = 'output2a <- FindCurve("CjaP", "EXT", c(output1$bifpoints[1,1:5], 10.0),  0.5, c(0.0, 1.0, 10.0, 50.0), parameters, c("par1", "9", "par2", "13", "EXTfun", "1", "EXTpar", "9", "noBP", "noLP", "silent"))'

eval(parse(text=cmd))

cmd = 'output2b <- FindCurve("CjaP", "EXT", c(output1$bifpoints[1,1:5], 10.0), -0.5, c(0.0, 1.0, 20.0, 50.0), parameters, c("par1", "9", "par2", "13", "EXTfun", "1", "EXTpar", "9", "noBP", "noLP", "silent"))'
eval(parse(text=cmd))

cat('\n\n\nContinuation of the the predator extinction boundary as a function of consumer mortality and interval\n\n');
cmd = 'output3a <- FindCurve("CjaP", "BP", c(output1$bifpoints[2,1:5], 10.0), 0.5, c(0.0, 1.0, 10.0, 50.0), parameters, c("par1", "9", "par2", "13", "silent"))'

eval(parse(text=cmd))

cmd = 'output3b <- FindCurve("CjaP", "BP", c(output1$bifpoints[2,1:5], 10.0), -0.5, c(0.0, 1.0, 10.0, 50.0), parameters, c("par1", "9", "par2", "13", "silent"))'
eval(parse(text=cmd))

parameters = c(0.1, 30.0, 1.0, 0.1, 0.5, 0.1, 0.015, 1.5, 3.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.32, 0.032, 0.5, 0.005, 0.0, 3.0, 0.0, 1.0, 0.0);

cat('\n\n\nTesting saddle-node bifurcation detection in the predator consumer-resource equilibrium\n\n');
cmd = 'output1 <- FindCurve("CjaP", "EQ", c(0.0, 27.0396, 0.590107, 0.132521, 1.84090), 1.0, c(0.0, 1.0), parameters, c("par1", "21", "par2", "13", "silent"))'

eval(parse(text=cmd))
cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

cat('\n\n\nContinuation of the saddle-node bifurcation as a function of predator mortality and interval\n\n');
cmd = 'output2 <- FindCurve("CjaP", "LP", output1$bifpoints[1,1:6], -0.2, c(0.0, 1.0, 20.0, 70.0), parameters, c("par1", "21", "par2", "13", "silent"))'

eval(parse(text=cmd))

cmd = 'output3 <- FindCurve("CjaP", "LP", output1$bifpoints[2,1:6], -0.2, c(0.0, 1.0, 20.0, 70.0), parameters, c("par1", "21", "par2", "13", "silent"))'
eval(parse(text=cmd))

cat('\n\n\nContinuation of the predator invasion boundary as a function of predator mortality and interval\n\n');
cmd = 'output4 <- FindCurve("CjaP", "BP", output1$bifpoints[3,1:6], -0.2, c(0.0, 1.0, 20.0, 70.0), parameters, c("par1", "21", "par2", "13", "silent"))'

eval(parse(text=cmd))

CleanFindCurve('F')
setwd(oldwd)

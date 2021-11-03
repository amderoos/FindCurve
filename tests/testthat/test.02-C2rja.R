devAskNewPage(ask = FALSE)
oldwd = getwd()

init = c(2.07278E-01, 2.07266E-01, 2.00000E+00, 1.00656E-05, 8.00686E-08);
parameters = c(1.0, 2.0, 2.0, 0.5, 0.1, 10.0, 1.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 10.0);

cat('\n\n\nConsumer-resource equilibrium as a function of R1max for default parameters\n\n');
cmd = 'output1 <- FindCurve("C2rja", "EQ", init, 0.1, c(0.0, 3.0), parameters, NULL, clean = TRUE, force = TRUE)'

eval(parse(text=cmd))
cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

CleanFindCurve('F')
setwd(oldwd)

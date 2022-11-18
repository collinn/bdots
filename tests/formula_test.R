# Need to formally test all this but not right now
# 
# ff1 <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1) + g3(m1, m2, m3, m4) + g4(r1)
# ff1 <- diffs(y, TrialType(M,W)) ~ Group(LI, TD)
# ff1 <- y ~ Group(LI, TD) + TrialType(M)
# 
# fit <- readRDS("../btest/eightgrpfit.rds")
# 
# ff1 <- y ~ Color(red, blue)
# l1 <- bootParser(ff1, fit)
# tt <- bootGroupSubset(l1, fit)
# 
# 
# 
# ff2 <- y ~ Color(red, blue) + Class(car)
# l2 <- bootParser(ff2, fit)
# tt2 <- bootGroupSubset(l2, fit)
# 
# ff3 <- diffs(y, Origin(foreign, domestic)) ~ Color(red, blue) + Class(car)
# l3 <- bootParser(ff3, fit)
# tt3 <- bootGroupSubset(l3, fit)
# 
# ff4 <- y ~ Class(car, truck) + Color(red) + Origin(foreign) 
# l4 <- bootParser(ff4, fit)
# tt4 <- bootGroupSubset(l4, fit)
# 
# 
# r1 <- bdotsBoot(ff1, fit)
# r2 <- bdotsBoot(ff2, fit)
# r3 <- bdotsBoot(ff3, fit)
# r4 <- bdotsBoot(ff4, fit)
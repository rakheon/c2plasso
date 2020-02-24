dat = read.csv("/Users/junkyungkim/Desktop/Actuarial/RMA/WorkersComp.csv")
dim(dat); head(dat)
dat1 = dat[seq(1,847,7),]
dim(dat1)
datre = reshape(dat[,-3], timevar = "YR", idvar = "CL", direction = "wide")
datre = as.matrix(datre)
x = t(datre[,-1])
x
x = x[,-c(18,22,65)]
S = cov(x)
S = cor(x)

library(glasso)
gl = glasso(S, rho = 0.7)
glwi = gl$wi
diag(glwi)=0
max(glwi); min(glwi)

plot(datre[4,-1])

# Test method on data used as example in Julie Josse's DenoiseR package
rm(list = ls())

library(devtools)
install_github("maryclare/LANOVA")
library(LANOVA)
library(RColorBrewer)
library(xtable)
library(leapp)

data(braintumors)

Y <- braintumors[["Y"]]
class <- braintumors[["class"]]

ett <- estTuneTest(Y)

# Reject null of no interactions
ett$test.stat > qnorm(0.95)
# P-value: want to compute Pr(z > test.stat)
exp(pnorm(ett$test.stat, log = TRUE, lower.tail = FALSE))

n <- nrow(Y); p <- ncol(Y)

H.1 <- diag(n) - rep(1, n)%*%t(rep(1, n))/n
H.2 <- diag(p) - rep(1, p)%*%t(rep(1, p))/p

fit <- estM(Y, ett)

m.f <- fit[["mu.f"]]
a.f <- fit[["a.f"]]
b.f <- fit[["b.f"]]
C.f <- fit[["C.f"]]
M.f <- fit[["M.f"]]

X <- matrix(NA, nrow = n*p, ncol = 1 + n + p - 2)
X[, 1] <- 1
X[, 2:p] <- (rep(1, nrow(Y))%x%diag(ncol(Y)))[, -1]
X[, (p + 1):ncol(X)] <- (diag(nrow(Y))%x%rep(1, ncol(Y)))[, -1]
XtX <- crossprod(X)
XtX.inv <- solve(XtX)
XXtX.inv <- tcrossprod(X, XtX.inv)
XXtX.invXt <- tcrossprod(XXtX.inv, X)

system.time(ipod.s <- IPOD(X = X, Y = c(Y), H = XXtX.invXt, method = "soft"))
c.ipod.s <- ipod.s$gamma
Theta.ipod.s <- matrix(lm(I(c(Y) - c.ipod.s)~X-1)$fitted.values + c.ipod.s, nrow = p, ncol = p)

ipod.h <- IPOD(X = X, Y = c(Y), H = H, method = "hard")
c.ipod.h <- ipod.h$gamma
Theta.ipod.h <- matrix(lm(I(c(Y) - c.ipod.h)~X-1)$fitted.values + c.ipod.h, nrow = p, ncol = p)



cols <- brewer.pal(11, "RdBu")
cols.tum <- c(cols[1:3], cols[length(cols) - 1])
cols <- ifelse(class == "A", cols.tum[1], ifelse(class == "O", cols.tum[2],
                                                 ifelse(class == "OA", cols.tum[3], cols.tum[4])))

Z.M <- t(H.1%*%M.f%*%H.2)
Z.Y <- t(H.1%*%Y%*%H.2)
colnames(Z.M) <- colnames(Z.Y) <- row.names(Y)
rownames(Z.M) <- rownames(Z.Y) <- colnames(Y)

par(mfrow = c(2, 1))
boxplot(Z.Y, col = cols)
boxplot(Z.M, col = cols)
# O3, GBM3, GBM4
par(mfrow = c(2, 1))
boxplot(t(Z.Y))
boxplot(t(Z.M), col = "pink")

zeros.C <- C.f != 0
zeros.C <- zeros.C[order(rowMeans(C.f != 0)), ]
class.new <- class[order(rowMeans(C.f != 0))]
cols.new <- ifelse(class.new == "A", cols.tum[1], 
                   ifelse(class.new == "O", cols.tum[2],
                          ifelse(class.new == "OA", cols.tum[3], cols.tum[4])))

plot(rowMeans(zeros.C), type = "n")
text(rowMeans(zeros.C), rownames(zeros.C), cex = 0.5, col = cols.new)
tail(rowMeans(zeros.C))

zeros.C <- C.f != 0
zeros.C <- zeros.C[, order(colMeans(C.f != 0))]

plot(colMeans(zeros.C), type = "n")
text(colMeans(zeros.C), colnames(zeros.C), cex = 0.5, col = cols.new)

vs <- names(tail(colMeans(zeros.C), 10))

par(mfrow = c(1, 2))
for (v in vs) {
  plot(Z.Y[v, ], type = "n", main = v)
  text(Z.Y[v, ], colnames(Z.Y), col = cols, cex = 0.5)
  
  plot(Z.M[v, ], type = "n", main = v)
  text(Z.M[v, ], colnames(Z.M), col = cols, cex = 0.5)
}


C.reo <- C.f[order(rowMeans(C.f != 0)), order(colMeans(C.f != 0))]
# Fix names that correspond to less common id's

# C.reo[C.reo == 0] <- NA
par(mfrow = c(1, 1))
# m <- max(abs(C.reo), na.rm = TRUE)
m <- c(quantile(abs(C.f[C.f!=0]), seq(0, 1, length.out = 6)))
breaks <- c(-m[length(m):1], m)
cols.scale <- brewer.pal(11, "RdBu")
cols.scale[6] <- "white"
cols.scale <- cols.scale[length(cols.scale):1]
q <- p - 50

pdf("josse_fig.pdf", family = "Times",
    width = 6.5, height = 4)
par(mfrow = c(1, 2))
par(cex.axis = 1)
par(cex.lab = 1)
par(cex.main = 1)
par(mar = c(3, 3, 3, 1))
par(mgp = c(1.75, 0.75, 0))
# par(oma = c(0.5, 0.5, 1, 0))

image(C.reo[1:n, 1:p], axes = FALSE,  breaks = breaks, col = cols.scale)
mtext("Gene", 2, line = 0)
mtext("Tumor", 1, line = 2)
axis(1, c(-0.1, 0:(nrow(C.reo)-1)/(nrow(C.reo)-1), 1.1), c("", rownames(C.reo), ""), las = 2, cex.axis = 0.45)
axis(2, c(0-0.01, 1.01), c("", ""))
# axis(3, c(-0.1, 0:(nrow(C.reo)-1)/(nrow(C.reo)-1), 1.1), c("", rownames(C.reo), ""), las = 2, cex.axis = 0.45)
eps <- 0.01
polygon(c(-0 - eps, 1+eps, 1+eps, 0-eps), c((q-1)/(p-1), (q-1)/(p-1), 1, 1))

cnames <- colnames(C.reo)
# cnames[cnames == "IMAGE.33267"] <- "IGFBP2"
image(C.reo[1:n, (p-50):p], axes = FALSE, breaks = breaks, col = cols.scale)
axis(1, c(-0.1, 0:(nrow(C.reo)-1)/(nrow(C.reo)-1), 1.1), c("", rownames(C.reo), ""), las = 2, cex.axis = 0.45)
axis(2, c(-0.1, 0:(p - q-1)/(p - q-1), 1.1), c("", cnames[(q + 1):p], ""), las = 2, cex.axis = 0.45)
# axis(3, c(-0.1, 0:(nrow(C.reo)-1)/(nrow(C.reo)-1), 1.1), c("", rownames(C.reo), ""), las = 2, cex.axis = 0.45)
mtext("Tumor", 1, line = 2)
dev.off()
# polygon(c(-0 - eps, 1+eps, 1+eps, 0-eps), c(0, 0, (p-q-1)/(p-1), (p-q-1)/(p-1)))

# Heterogeneity among glioblastoma?, e.g. DLLP

# This gives the bottom left corner (to check indices)
t(C.reo)[p:(p-5), (n-5):n]

# de Tayrac etal (these are already ordered by sparsity)
tyr.genes <- c("MSN", "RUNX1", "PPP3CB", "ARPP.19", "H08563", "C9orf48",
               "RTN3", "CLIC1", "WASF1", "S100A11", "UBA52", "VAMP2", 
               "AA398420", "RALY", "X37864", "PDXP", "EMP3", 
               "AA281932", "ASPA", "PDPN")

par(mfrow = c(1, 1))
image(C.reo[, tyr.genes][1:n, 1:length(tyr.genes)], axes = FALSE, breaks = breaks, col = cols.scale)
axis(1, c(-0.1, 0:(nrow(C.reo)-1)/(nrow(C.reo)-1), 1.1), c("", rownames(C.reo), ""), las = 2, cex.axis = 0.45)
# axis(2, c(-0.1, 0:(length(tyr.genes) - 1)/(length(tyr.genes) - 1), 1.1), c("", colnames(C.reo)[colnames(C.reo) %in% tyr.genes], ""), las = 2, cex.axis = 0.45)
axis(2, c(-0.1, 0:(length(tyr.genes) - 1)/(length(tyr.genes) - 1), 1.1), c("", tyr.genes, ""), las = 2, cex.axis = 0.45)

par(mfrow = c(2, 2))
for (tyg in c("ASPA", "PDPN")) {
  
  plot(Z.Y[tyg, ], type = "n", main = tyg)
  text(Z.Y[tyg, ], colnames(Z.Y), col = ifelse(class == "GBM", "blue", "red"),
       cex = 0.5)
  
  plot(Z.M[tyg, ], type = "n", main = tyg)
  text(Z.M[tyg, ], colnames(Z.M), col = ifelse(class == "GBM", "blue", "red"),
       cex = 0.5)
}
for (tyg in c("MSN", "RUNX1")) {
  
  plot(Z.Y[tyg, ], type = "n", main = tyg)
  text(Z.Y[tyg, ], colnames(Z.Y), col = ifelse(class == "GBM", "blue", "red"),
       cex = 0.5)
  
  plot(Z.M[tyg, ], type = "n", main = tyg)
  text(Z.M[tyg, ], colnames(Z.M), col = ifelse(class == "GBM", "blue", "red"),
       cex = 0.5)
}

par(mfrow = c(1, 2))
for (mcg in cnames[length(cnames):(length(cnames) - 20)]) {
  
  plot(Z.Y[mcg, ], type = "n", main = mcg)
  text(Z.Y[mcg, ], colnames(Z.Y),
       cex = 0.5, col = cols)
  
  plot(Z.M[mcg, ], type = "n", main = mcg)
  text(Z.M[mcg, ], colnames(Z.M), col = cols,
       cex = 0.5)
}

three <- t(C.reo[(nrow(C.reo) - 2):nrow(C.reo), apply(C.reo[(nrow(C.reo) - 2):nrow(C.reo), ], 2, function(x) {min(abs(x)) > 0})])
# Try to figure out what's going on with the top 3 glioblastomas
top.genes <- row.names(three) # All of these perfectly separate top 3 from rest
top.genes <- c("FCGR2B", # Immunological risk factor? Overexpressed when prog worse 
               "CACNB3", # Couldn't find anything about this
               "RTN3", # Differentially expressed in different glioblastoma subtypes
               "F13A1", 
               "NY.SAR.48",
               "HMOX1", # overexpression -> poor prognosis, Ghosh2016a
               "AA281932", # ???, found description "pre-B-cell colony-enhancing factor" in Ross and Perou
               "PLAUR", 
               "TGFBI", # Overexpressed in primary glioblastoma
               "AI262682", # AKA FN1: http://www2.stat.duke.edu/~sayan/ashley/Vascular/Inflammatory_Response_Pathway.html, overexpressed in primary
               "PEG3", # Downreg in ovarian cancer
               "S100A11") # Overexpressed in primary glio
for (top in top.genes) {
  plot(Y[, top], type = "n", main = top)
  text(Y[,  top], row.names(Y), col = ifelse(row.names(Y) %in% c("GBM3", "GBM4", "GBM30"), "red", "black"))
  
}

C.f[, "IL10"]

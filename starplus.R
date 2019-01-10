rm(list = ls())

library(RColorBrewer)
library(devtools)
install_github("maryclare/LANOVA")
library(LANOVA)

data(fMRI)

Y <- fMRI[["Y"]]
map <- fMRI[["map"]]
xs <- map[, 1]; ys <- map[, 2]; zs <- map[, 3]

# Obtain tuning parameter estimates
ett <- estTuneTest(Y, lowpen = TRUE)
# Reject null of no interactions
ett$test.stat > qnorm(0.95)
# P-value: want to compute Pr(z > test.stat)
exp(pnorm(ett$test.stat, log = TRUE, lower.tail = FALSE))


fit <- estM(Y, ett)

n<-dim(Y)[1] ; p<-dim(Y)[2] ; q <- dim(Y)[3]

m.1 <- t(rep(1, n)/n)
m.2 <- t(rep(1, p)/p)
m.3 <- t(rep(1, q)/q)

H.1 <- diag(n) - rep(1, n)%*%t(rep(1, n))/n
H.2 <- diag(p) - rep(1, p)%*%t(rep(1, p))/p
H.3 <- diag(q) - rep(1, q)%*%t(rep(1, q))/q

mu.f <- fit$mu.f
a.f <- fit$a.f
b.f <- fit$b.f
d.f <- fit$d.f
E.f <- fit$E.f
F.f <- fit$F.f
G.f <- fit$G.f
C.f <- fit$C.f
M.f <- fit$M.f

mu.l <- atrans(fit$M.f, list(m.1, m.2, m.3))[, , 1]
a.l <- atrans(fit$M.f, list(H.1, m.2, m.3))[, , 1]
b.l <- atrans(fit$M.f, list(m.1, H.2, m.3))[, , 1]
d.l <- atrans(fit$M.f, list(m.1, m.2, H.3))[1, , ]
E.l <- atrans(fit$M.f, list(H.1, H.2, m.3))[, , 1]
F.l <- atrans(fit$M.f, list(H.1, m.2, H.3))[, 1, ]
G.l <- atrans(fit$M.f, list(m.1, H.2, H.3))[1, , ]
C.l <- atrans(fit$M.f, list(H.1, H.2, H.3))

mu.o <- atrans(Y, list(m.1, m.2, m.3))[, , 1]
a.o <- atrans(Y, list(H.1, m.2, m.3))[, , 1]
b.o <- atrans(Y, list(m.1, H.2, m.3))[, , 1]
d.o <- atrans(Y, list(m.1, m.2, H.3))[1, , ]
E.o <- atrans(Y, list(H.1, H.2, m.3))[, , 1]
F.o <- atrans(Y, list(H.1, m.2, H.3))[, 1, ]
G.o <- atrans(Y, list(m.1, H.2, H.3))[1, , ]
C.o <- atrans(Y, list(H.1, H.2, H.3))

sum(c(mu.f, a.f, b.f, d.f,
      c(E.f), c(F.f), c(G.f), c(C.f)) != 0)

sum(c(mu.f, a.f, b.f, d.f,
      c(E.f), c(F.f), c(G.f), c(C.f)) != 0)/prod(dim(Y))

mean(a.f == 0)
mean(b.f == 0)
mean(d.f == 0)

mean(E.f == 0)
mean(F.f == 0)
mean(G.f == 0)

mean(C.f == 0)

plot(c(F.o), c(F.l))

cols <- brewer.pal(9, "Reds")
breaks <- seq(0, 1, length.out = 10)

F.f.mean <- apply(F.f != 0, 2, mean)

F.cols.for.plot <- ifelse(F.f.mean <= breaks[2], cols[1],
                        ifelse(F.f.mean <= breaks[3], cols[2],
                               ifelse(F.f.mean <= breaks[4], cols[3],
                                      ifelse(F.f.mean <= breaks[5], cols[4],
                                             ifelse(F.f.mean <= breaks[6], cols[5],
                                                    ifelse(F.f.mean <= breaks[7], cols[6],
                                                           ifelse(F.f.mean <= breaks[8], cols[7],
                                                                  ifelse(F.f.mean <= breaks[9], cols[8], cols[9]))))))))

C.f.mean <- apply(C.f != 0, 3, mean)

C.cols.for.plot <- ifelse(C.f.mean <= breaks[2], cols[1],
                        ifelse(C.f.mean <= breaks[3], cols[2],
                               ifelse(C.f.mean <= breaks[4], cols[3],
                                      ifelse(C.f.mean <= breaks[5], cols[4],
                                             ifelse(C.f.mean <= breaks[6], cols[5],
                                                    ifelse(C.f.mean <= breaks[7], cols[6],
                                                           ifelse(C.f.mean <= breaks[8], cols[7],
                                                                  ifelse(C.f.mean <= breaks[9], cols[8], cols[9]))))))))


pdf("fmri.pdf", family = "Times",
    width = 6, height = 2)
par(cex.axis = 1)
par(cex.lab = 1)
par(cex.main = 1)
par(mgp = c(1.75, 0.75, 0))
par(mfrow = c(2, 8))
par(mar = rep(0.25, 4))
par(oma = c(0, 2, 0, 0))
for (i in 1:8) {
  plot(-xs[zs == i], -ys[zs == i], col = F.cols.for.plot,
       pch = 15, xlab = "", ylab = "", axes = FALSE, cex = 0.5)
}
mtext(expression(paste("                           ", hat(F), sep = "")), 2, outer = TRUE, line = 0.5)
for (i in 1:8) {
  plot(-xs[zs == i], -ys[zs == i], col = C.cols.for.plot,
       pch = 15, xlab = "", ylab = "", axes = FALSE, cex = 0.5)
}
mtext(expression(paste(hat(C), "                            ", sep = "")), 2, outer = TRUE, line = 0.5)
dev.off()


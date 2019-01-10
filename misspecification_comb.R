rm(list = ls())

library(RColorBrewer)
cols <- brewer.pal(11, "RdBu")
cols[6] <- cols[3]
cols <- cols[c(2, 6, length(cols) - 1)]

k <- 2

pdf(paste("misspecification_comb_", k, ".pdf", sep = ""), family = "Times",
    width = 6, height = 4)
cex.lab <- 0.7
par(mfrow = c(2, 4))
par(cex.axis = 1)
par(cex.lab = 1)
par(cex.main = 1)
par(mar = c(3, 1, 1, 1))
par(mgp = c(1.75, 0.75, 0))
par(oma = c(0.5, 2.5, 1, 0))

load(paste("misspecification_gnorm_", k, ".RData", sep = ""))
sig.sq.ma <- sig.sq.zs <- 
  mse.uni <- mse.sur <- mse.lan <- mse.mle <- mse.add <- mse.low.1 <- mse.low.5 <- mse.ipod.h <- mse.ipod.s <- 
  mae.uni <- mae.sur <- mae.lan <- mae.mle <- mae.add <- mae.low.1 <- mae.low.5 <- mae.ipod.h <- mae.ipod.s <- 
  mzo.uni <- mzo.sur <- mzo.lan <- mzo.mle <- mzo.add <- mzo.low.1 <- mzo.low.5 <- mzo.ipod.h <- mzo.ipod.s <- 
  lanova.time <- ipod.s.time <- ipod.h.time <- array(dim = c(length(s2cs), sims, length(prs)))

for (i in 1:(length(prs))) {
  mse.lan[, , i] <- sim.prs[[i]]$mse.lan
  mse.mle[, , i] <- sim.prs[[i]]$mse.mle
  mse.uni[, , i] <- sim.prs[[i]]$mse.uni
  mse.sur[, , i] <- sim.prs[[i]]$mse.sur
  mse.add[, , i] <- sim.prs[[i]]$mse.add
  mse.low.1[, , i] <- sim.prs[[i]]$mse.low.1
  mse.low.5[, , i] <- sim.prs[[i]]$mse.low.5
  mse.ipod.s[, , i] <- sim.prs[[i]]$mse.ipod.s
  mse.ipod.h[, , i] <- sim.prs[[i]]$mse.ipod.h
  
  sig.sq.ma[, , i] <- sim.prs[[i]]$sig.sq.ma
  sig.sq.zs[, , i] <- sim.prs[[i]]$sig.sq.zs
}

ylim <- c(-1.5, 1.5)
xlim <- range(prs)
plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.mle[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,  main = "Least Squares")
mtext("Log Relative Risk", 2, line = 2, cex = cex.lab)
mtext(expression(q[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.mle[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.mle[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.mle[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.mle[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.mle[3, , ])), col = cols[3])

abline(h = 0, lty = 3)

points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.add[1, , ])), col = cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.add[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.add[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.add[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.add[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.add[3, , ])), col = cols[3], lty = 2)

legend("topleft", legend = c("MLE", "Add."),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")

plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.1[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,
     main = "Low Rank Int.")
mtext("", 2, line = 2, cex = cex.lab)
mtext(expression(q[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.1[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.1[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.1[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.1[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.1[3, , ])), col = cols[3])
points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.5[1, , ])), col= cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.5[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.5[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.5[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.5[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.5[3, , ])), col = cols[3], lty = 2)
abline(h = 0, lty = 3)

legend("topright", legend = c("Rank 1", "Rank 5"),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")

plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.s[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,
     main = "IPOD")
mtext("", 2, line = 2, cex = cex.lab)
mtext(expression(q[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.s[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.s[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.s[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.s[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.s[3, , ])), col = cols[3])
points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.h[1, , ])), col= cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.h[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.h[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.h[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.h[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.h[3, , ])), col = cols[3], lty = 2)
abline(h = 0, lty = 3)

legend("topright", legend = c("Soft Thresh.", "Hard Thresh."),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")


plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.uni[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,
     main = "Minimax Int.")
mtext("", 2, line = 2, cex = cex.lab)
mtext(expression(q[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.uni[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.uni[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.uni[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.uni[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.uni[3, , ])), col = cols[3])

points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.sur[1, , ])), col= cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.sur[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.sur[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.sur[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.sur[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.sur[3, , ])), col = cols[3], lty = 2)

abline(h = 0, lty = 3)

legend("topright", legend = c("Universal", "SURE"),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")

legend.text <- c(expression(paste(sigma[c]^2, "=", sigma[z]^2, "/2", sep = "")),
                 expression(paste(sigma[c]^2, "=", sigma[z]^2, sep = "")),
                 expression(paste(sigma[c]^2, "=2", sigma[z]^2, "", sep = "")))

legend("bottomright", legend = legend.text,
       pch = c(16, 17, 18),
       col = c(cols[1], cols[2], cols[3]), bg = "white", cex = 1,
       bty = "n")

load(paste("misspecification_", k, ".RData", sep = ""))
sig.sq.ma <- sig.sq.zs <- 
  mse.uni <- mse.sur <- mse.lan <- mse.mle <- mse.add <- mse.low.1 <- mse.low.5 <- mse.ipod.h <- mse.ipod.s <- 
  mae.uni <- mae.sur <- mae.lan <- mae.mle <- mae.add <- mae.low.1 <- mae.low.5 <- mae.ipod.h <- mae.ipod.s <- 
  mzo.uni <- mzo.sur <- mzo.lan <- mzo.mle <- mzo.add <- mzo.low.1 <- mzo.low.5 <- mzo.ipod.h <- mzo.ipod.s <- 
  lanova.time <- ipod.s.time <- ipod.h.time <- array(dim = c(length(s2cs), sims, length(prs)))

for (i in 1:(length(prs))) {
  mse.lan[, , i] <- sim.prs[[i]]$mse.lan
  mse.mle[, , i] <- sim.prs[[i]]$mse.mle
  mse.uni[, , i] <- sim.prs[[i]]$mse.uni
  mse.sur[, , i] <- sim.prs[[i]]$mse.sur
  mse.add[, , i] <- sim.prs[[i]]$mse.add
  mse.low.1[, , i] <- sim.prs[[i]]$mse.low.1
  mse.low.5[, , i] <- sim.prs[[i]]$mse.low.5
  mse.ipod.s[, , i] <- sim.prs[[i]]$mse.ipod.s
  mse.ipod.h[, , i] <- sim.prs[[i]]$mse.ipod.h
  
  sig.sq.ma[, , i] <- sim.prs[[i]]$sig.sq.ma
  sig.sq.zs[, , i] <- sim.prs[[i]]$sig.sq.zs
}

xlim <- range(prs)
plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.mle[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,  main = "")
mtext("Log Relative Risk", 2, line = 2, cex = cex.lab)
mtext(expression(pi[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.mle[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.mle[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.mle[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.mle[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.mle[3, , ])), col = cols[3])

abline(h = 0, lty = 3)

points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.add[1, , ])), col = cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.add[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.add[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.add[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.add[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.add[3, , ])), col = cols[3], lty = 2)

legend("topleft", legend = c("MLE", "Add."),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")

plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.1[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,
     main = "")
mtext("", 2, line = 2, cex = cex.lab)
mtext(expression(pi[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.1[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.1[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.1[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.1[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.1[3, , ])), col = cols[3])
points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.5[1, , ])), col= cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.low.5[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.5[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.low.5[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.5[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.low.5[3, , ])), col = cols[3], lty = 2)
abline(h = 0, lty = 3)

legend("topright", legend = c("Rank 1", "Rank 5"),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")

plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.s[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,
     main = "")
mtext("", 2, line = 2, cex = cex.lab)
mtext(expression(pi[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.s[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.s[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.s[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.s[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.s[3, , ])), col = cols[3])
points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.h[1, , ])), col= cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.ipod.h[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.h[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.ipod.h[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.h[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.ipod.h[3, , ])), col = cols[3], lty = 2)
abline(h = 0, lty = 3)

legend("topright", legend = c("Soft Thresh.", "Hard Thresh."),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")


plot(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.uni[1, , ])), xlab = "", ylab = "",
     ylim = ylim, col = cols[1], pch = 16, xlim = xlim,
     main = "")
mtext("", 2, line = 2, cex = cex.lab)
mtext(expression(pi[c]), 1, line = 2, cex = cex.lab)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.uni[1, , ])), col= cols[1])
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.uni[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.uni[2, , ])), col= cols[2])
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.uni[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.uni[3, , ])), col = cols[3])

points(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.sur[1, , ])), col= cols[1], pch = 16)
lines(prs, log(colMeans(mse.lan[1, , ])/colMeans(mse.sur[1, , ])), col= cols[1], lty = 2)
points(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.sur[2, , ])), col= cols[2], pch = 17)
lines(prs, log(colMeans(mse.lan[2, , ])/colMeans(mse.sur[2, , ])), col= cols[2], lty = 2)
points(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.sur[3, , ])), col = cols[3], pch = 18)
lines(prs, log(colMeans(mse.lan[3, , ])/colMeans(mse.sur[3, , ])), col = cols[3], lty = 2)

abline(h = 0, lty = 3)

legend("topright", legend = c("Universal", "SURE"),
       lty = c(1, 2), bg = "white", cex = 1,
       bty = "n")

legend.text <- c(expression(paste(tau[c]^2, "=", sigma[z]^2, "/2", sep = "")),
                 expression(paste(tau[c]^2, "=", sigma[z]^2, sep = "")),
                 expression(paste(tau[c]^2, "=2", sigma[z]^2, "", sep = "")))
legend("bottomright", legend = legend.text,
       pch = c(16, 17, 18),
       col = c(cols[1], cols[2], cols[3]), bg = "white", cex = 1,
       bty = "n")

dev.off()

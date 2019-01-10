rm(list = ls())

library(RColorBrewer)

cols.ss <- brewer.pal(11, "RdBu")
cols.ss[6] <- "black"

rs.ss <- seq(0, 2, length.out = 11)
prs <- seq(0, 1, 0.01)

analytic.rates.ss <- array(dim = c(length(rs.ss), length(prs), 2)) 

ns.ss <- c(100, 1000)
for (r in rs.ss) {
  for (pr in prs) {
    
    for (n in ns.ss) {
    mu.test <- sqrt(3*n/8)*pr*(1 - pr)*r^2/(pr*r + 1)^2
    var.test <- 1 + pr*(1 - pr)*((20*pr^2 - 28*pr + 35)*r^4 + 16*(5 - pr)*r^3 + 72*r^2)/(8*(pr*r + 1)^4)
    analytic.rates.ss[which(r == rs.ss), which(pr == prs), which(n == ns.ss)] <- 1 - pnorm((qnorm(0.95)-mu.test)/sqrt(var.test))
    }
  }
}

cols.lap <- brewer.pal(9, "Blues")
cols.lap <- c(cols.lap[2:9], "black")

rs.lap <- seq(0, 2, by = 0.01)
ns.lap <- (1:9)*100 + 100
analytic.rates.lap <- matrix(nrow = length(rs.lap), ncol = length(ns.lap))

for (r in rs.lap) {
  for (n in ns.lap) {
    
    mu.test <- sqrt(n/24)*3*r^2/(1 + r)^2
    
    var.test <- 1 + (68*r^4 + 38*r^3 + 9*r^2)/(1 + r)^4
    
    analytic.rates.lap[which(r == rs.lap), which(n == ns.lap)] <- 1 - pnorm((qnorm(0.95)-mu.test)/sqrt(var.test))
    
  }
}

pdf("power.pdf", family = "Times",
    width = 6, height = 2.5)
par(cex.axis = 1)
par(cex.lab = 1)
par(cex.main = 1)
par(mfrow = c(1, 3))
# par(mar = c(5.1,4.1,4.1,2.1) - c(0, 0.5, 0, 2))

plot(rs.lap, analytic.rates.lap[, 1], type = "n", ylim = range(analytic.rates.lap),
     xlab = "", ylab = "",
     axes = FALSE)

for (i in 1:length(ns.lap)) {
  
  lines(rs.lap, analytic.rates.lap[, i], lty = 1, col = cols.lap[i])
}
mtext("Power", 2 , line = 2, cex = 0.75)
mtext(expression(paste(phi^2, "=", sigma[c]^2/sigma[z]^2, sep = "")), 1 , line = 2.5, cex = 0.75)
mtext("Laplace", 3 , line = 1, cex = 0.75)
axis(2, pretty(as.vector(analytic.rates.ss)))
axis(1, c(-1, pretty(rs.lap)), labels = c("", pretty(rs.lap)))
legend("bottomright", lty = rep(1, 2), legend = ns.lap[c(1, length(ns.lap))], col = cols.lap[c(1, length(cols.lap))], title = "Values of np", bty = "n",
       cex = 1, bg = "white")

for (j in 1:2) {
plot(prs, analytic.rates.ss[1, , j], type = "n", ylim = range(analytic.rates.ss),
     xlab = "", ylab = "", axes = FALSE, main = )
  mtext(expression(pi[c]), 1 , line = 2.5, cex = 0.75)
  mtext(paste("Spike-and-slab, np=", ns.ss[j], sep = ""), 3 , line = 1, cex = 0.75)
for (i in 1:length(rs.ss)) {
  lines(prs, analytic.rates.ss[i, , j], lty = 1, col = cols.ss[i])
}
  axis(2, pretty(as.vector(analytic.rates.ss)), labels = rep("", 6))
  if (j == 1) {
    mid <- median(1:length(cols.ss))
    legend("topright", lty = rep(1, 3), col = cols.ss[c(1, mid, length(cols.ss))], legend = rs.ss[c(1, mid, length(rs.ss))],
           title = expression(paste("Values of ", phi^2, "=", tau[c]^2/sigma[z]^2, sep = "")), bty = "n", cex = 1)
    mtext("", side = 2, line = 2, cex = 1)
  }
  axis(1, c(-1, pretty(prs)), labels = c("", pretty(prs)))
}
dev.off()

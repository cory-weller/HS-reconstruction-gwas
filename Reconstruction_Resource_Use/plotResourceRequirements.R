library(ggplot2)
library(data.table)
library(ggthemes)

dat <- fread('RABBIT_resource_requirements.tab')
dat.averaged <- dat[, list(runtimeMin = mean(runtimeMin), memoryGB = mean(memoryGB)), by=list(nSNPs, nFounders)]


# memory usage

fit.memory <- summary(lm(data=dat.averaged, memoryGB ~ I(nFounders^4):nSNPs))

#   Call:
#   lm(formula = memoryGB ~ I(nFounders^4):nSNPs, data = dat.averaged)
#
#   Residuals:
#        Min       1Q   Median       3Q      Max
#   -0.28973 -0.03775  0.03622  0.04675  0.07540
#
#   Coefficients:
#                         Estimate Std. Error  t value Pr(>|t|)
#   (Intercept)          3.196e-02  9.117e-03    3.505  0.00072 ***
#   I(nFounders^4):nSNPs 7.867e-09  4.185e-12 1879.921  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#   Residual standard error: 0.07221 on 88 degrees of freedom
#   Multiple R-squared:      1,     Adjusted R-squared:      1
#   F-statistic: 3.534e+06 on 1 and 88 DF,  p-value: < 2.2e-16

# Estimate runtime with model coefficients

dat.averaged[, memoryGBest := nSNPs*(nFounders^4)*7.867e-09  + 3.196e-02]

# format data in long form for plotting
dat.averaged.mem.long <- melt(dat.averaged, measure.vars=c("memoryGB", "memoryGBest"))

memoryModel5000 <- function(x) {
    S <- 5000
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel4500 <- function(x) {
    S <- 4500
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel4000 <- function(x) {
    S <- 4000
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel3500 <- function(x) {
    S <- 3500
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel3000 <- function(x) {
    S <- 3000
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel2500 <- function(x) {
    S <- 2500
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel2000 <- function(x) {
    S <- 2000
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel1500 <- function(x) {
    S <- 1500
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel1000 <- function(x) {
    S <- 1000
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}
memoryModel500 <- function(x) {
    S <- 500
    return(S*(x^4)*7.867e-09  + 3.196e-02)
}


runtimeModel5000 <- function(x) {
    S <- 5000
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel4500 <- function(x) {
    S <- 4500
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel4000 <- function(x) {
    S <- 4000
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel3500 <- function(x) {
    S <- 3500
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel3000 <- function(x) {
    S <- 3000
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel2500 <- function(x) {
    S <- 2500
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel2000 <- function(x) {
    S <- 2000
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel1500 <- function(x) {
    S <- 1500
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel1000 <- function(x) {
    S <- 1000
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}
runtimeModel500 <- function(x) {
    S <- 500
    return((1.189e-03*x^2 + 1.038e-06*S*x^2 + 2.649e-04*S )^(2))
}


# Plot memory usage
#  + geom_line(data=dat.averaged.mem.long[variable=="memoryGBest"], aes(x=nFounders, y=value, group=nSNPs), linetype="solid", alpha=0.3) +

g.memory <- ggplot() + geom_point(data=dat.averaged.mem.long[variable=="memoryGB"], aes(x=nFounders, y=value)) +
theme_few(12) +
scale_x_continuous(breaks=c(4,8,12,16,20,24,28,32,36), limits=c(4,42)) +
geom_text(dat.averaged.mem.long[nFounders==36 & variable=="memoryGBest"], mapping=aes(x=40, y=value+1, label=paste(nSNPs, "SNPs", sep=" ")), size=2.5) +
geom_text(data.table(x=18, y=60, label=("GB=7.867×10⁻⁹(SN⁴) + 0.0316\nR² = 1")), mapping=aes(x=x, y=y, label=label), size=2.5) +
labs(x="Number of Founding Haplotypes", y="Peak Memory Usage (GB)") +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel5000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel4500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel4000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel3500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel3000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel2500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel2000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel1500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel1000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=memoryModel500) +
ylim(0,75)

# end memory



ggsave(g.memory, file="memory.model.png", width=12, height=12, units="cm")
ggsave(g.memory, file="memory.model.svg", width=12, height=12, units="cm")


# Runtime
fit.runtime <- summary(lm(data=dat.averaged[nFounders>8], I(runtimeMin^(0.5)) ~ I(nFounders^2) * nSNPs + 0))    # R-squared = 0.9995

#   Call:
#   lm(formula = I(runtimeMin^(0.5)) ~ I(nFounders^2) * nSNPs + 0,
#       data = dat.averaged[nFounders > 8])
#
#   Residuals:
#        Min       1Q   Median       3Q      Max
#   -0.25060 -0.04379 -0.00500  0.04211  0.18058
#
#   Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
#   I(nFounders^2)       1.189e-03  3.166e-05   37.55   <2e-16 ***
#   nSNPs                2.649e-04  6.820e-06   38.84   <2e-16 ***
#   I(nFounders^2):nSNPs 1.038e-06  1.284e-08   80.86   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#   Residual standard error: 0.09177 on 67 degrees of freedom
#   Multiple R-squared:  0.9995,    Adjusted R-squared:  0.9995
#   F-statistic: 4.316e+04 on 3 and 67 DF,  p-value: < 2.2e-16
#

# Estimate runtime with model coefficients
dat.averaged[, runtimeMinest := (1.189e-03*nFounders^2 + 1.038e-06*nSNPs*nFounders^2 + 2.649e-04*nSNPs )^(2) ]

# format data in long form for plotting
dat.averaged.runtime.long <- melt(dat.averaged, measure.vars=c("runtimeMin", "runtimeMinest"))

# Plot runtime
# geom_line(data=dat.averaged.runtime.long[variable=="runtimeMinest"], aes(x=nFounders, y=value, group=nSNPs), linetype="solid", alpha=0.3) +

g.runtime <- ggplot() + geom_point(data=dat.averaged.runtime.long[variable=="runtimeMin"], aes(x=nFounders, y=value)) +
theme_few(12) +
scale_x_continuous(breaks=c(4,8,12,16,20,24,28,32,36), limits=c(4,42)) +
geom_text(dat.averaged.runtime.long[nFounders==36 & variable=="runtimeMin"], mapping=aes(x=40, y=value+1, label=paste(nSNPs, "SNPs", sep=" ")), size=2.5) +
geom_text(data.table(x=18, y=85, label=("Minutes=1.038×10⁻⁶(SN²) + 1.189×10⁻³(N²) + 2.649×10⁻⁴(S)\nR² = 0.9995")), mapping=aes(x=x, y=y, label=label), size=2.5) +
labs(x="Number of Founding Haplotypes", y="Runtime (Minutes)") +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel5000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel4500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel4000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel3500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel3000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel2500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel2000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel1500) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel1000) +
stat_function(data=data.table(x=4:36), aes(x), alpha=0.3, fun=runtimeModel500) +
ylim(0,95)

ggsave(g.runtime, file="runtime.model.svg", width=12, height=12, units="cm")

# end runtime

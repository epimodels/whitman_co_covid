### Transmission Parameter Estimation - Two-population SEIR Model ###


## ----model-construct----------------------------------------------------------
library(foreach)
library(iterators)
library(parallel)
library(rngtools)
library(doParallel)

library(doRNG)
registerDoParallel()
registerDoRNG(2488820)

library(pomp)
library(ggpubr) # use ggarrange() to combine separate figures in a single plot.
library(ggrepel) # use geom_text_repel() to place labels that don't overlap each other
library(bayestestR) # runs the ci function to get credible intervals
library(coda)
library(MCMCglmm) # processes the mcmc lists [requires Matrix, coda, ape packages]
library(MASS)
library(reshape2)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)
library(tibble)
library(readr)
library(stringr)
library(forcats)

set.seed(420)

# dir.create("tmp")
options(pomp_cdir="./tmp")


############# State and Parameter Definitions ############# 

# S_p notation should now be S_c
# S_w notation should now be S_u

# S_c - community susceptibles
# S_u - university susceptibles
# E_c - community exposed
# E_u - university exposed
# I_c - community infecteds
# I_u - university infecteds
# R_c - community recovered
# R_u - university recovered

### Accumulator variables of those who could be tested, not those infected:
# Z_c - community  
# Z_u - university


### Transmission Parameters

# beta_c  --> community transmission

# Beta_u  --> university transmission

# beta_m  --> community-university transmission

### Additional parameters

# alpha  --> time to infectiousness
# gamma  --> recovery rate of infecteds
# phi    --> testing probability (phi*accumulator variables) 
              # phi tests a proportion of those who are infected and eligible for testing
# k      --> variance constant


############# Csnippet SEIR Model ############# 

SEIR.compart <- Csnippet("
  double dN_S_cE_c = rbinom(S_c,1-exp((-beta_c/100000*I_c-beta_m/100000*I_u)*dt));
  double dN_S_uE_u = rbinom(S_u,1-exp((-beta_u/100000*I_u-beta_m/100000*I_c)*dt));

  double dN_E_cI_c = rbinom(E_c,1-exp(-alpha*dt));
  double dN_E_uI_u = rbinom(E_u,1-exp(-alpha*dt));
  
  double dN_I_cR_c = rbinom(I_c,1-exp(-gamma*dt));
  double dN_I_uR_u = rbinom(I_u,1-exp(-gamma*dt));
  
  S_c -= dN_S_cE_c;
  S_u -= dN_S_uE_u;
  
  E_c += dN_S_cE_c - dN_E_cI_c;
  E_u += dN_S_uE_u - dN_E_uI_u;
  
  I_c += dN_E_cI_c - dN_I_cR_c;
  I_u += dN_E_uI_u - dN_I_uR_u;
  
  R_c += dN_I_cR_c;
  R_u += dN_I_uR_u;
  
  Z_c += dN_E_cI_c;
  Z_u += dN_E_uI_u;
")

SEIR_rinit <- Csnippet("
  S_c = N_c-E_c0-I_c0-R_c0;
  S_u = N_u-E_u0-I_u0-R_u0;
  E_c = E_c0;
  E_u = E_u0;
  I_c = I_c0;
  I_u = I_u0;
  R_c = R_c0;
  R_u = R_u0;
  Z_c = 0;
  Z_u = 0;
")


dmeas <- Csnippet("if (beta_c<0||beta_u<0||beta_m<0||phi<0||k<0) {
  lik = (give_log) ? R_NegInf : 0.0;
} else {
  lik=dnbinom_mu(reports_c,k,(phi*Z_c),give_log)+dnbinom_mu(reports_u,k,(phi*Z_u),give_log);
}")

#dmeas gives the log likelihood, give_log default == 1


rmeas <- Csnippet("
  reports_c = rnbinom_mu(k,phi*Z_c);
  reports_u = rnbinom_mu(k,phi*Z_u);
  ")


############# Simulate Data #############

Week <- seq(31, 53, 1)
Week <- as.data.frame(Week)
Week$reports_c <- rep(NA, length(Week))
Week$reports_u <- rep(NA, length(Week$reports_c))

accumvars=c("Z_c","Z_u")
statenames=c("S_c","S_u","E_c","E_u","I_c","I_u","R_c","R_u","Z_c","Z_u")
paramnames=c("N_c","N_u",
             "E_c0","I_c0","R_c0",
             "E_u0","I_u0","R_u0",
             "beta_c","beta_u","beta_m",
             "alpha","gamma","phi","k")
true.theta=c(N_c=20785,N_u=14254,
             E_c0=0,I_c0=124,R_c0=19,
             E_u0=0,I_u0=81,R_u0=15,
             alpha=1/3.59,gamma=1/5,
             # Unknown Parameters:
             beta_c=2.1,
             beta_u=26,
             beta_m=0.04,
             phi=0.12,
             k=2)

Week %>%
  pomp(
    times="Week",t0=30,
    rprocess=euler(SEIR.compart,delta.t=1/7),
    rinit=SEIR_rinit,
    rmeasure=rmeas,              # cannot simulate without rmeas
    dmeasure=dmeas,              # only need dmeas for the pfilter
    accumvars=accumvars,
    statenames=statenames,
    paramnames=paramnames
  ) -> SEIR

SEIR %>%
  simulate(params=true.theta,
           nsim=1,
           format="data.frame",
           include.data=FALSE
  ) -> SEIR_sims

#View(SEIR_sims)

### Plot the simulated data
sim.data <- ggplot(SEIR_sims, aes(x=Week)) + theme_minimal() + ylab("Case Reports") + xlab("Week") +
  geom_line(aes(y=reports_c), color="darkred") +
  geom_line(aes(y=reports_u), color="steelblue") +
  guides(color="none") + ggtitle("Measurement Model")

plot(sim.data)

numeric.sol <- ggplot(SEIR_sims, aes(x=Week, group=.id)) + theme_minimal() + ylab("N") + xlab("Week") +
  geom_line(aes(y = S_c), color = "darkred", linetype=2) + 
  geom_line(aes(y = S_u), color = "darkred") + 
  geom_line(aes(y = E_c), color = "steelblue", linetype=2) + 
  geom_line(aes(y = E_u), color = "steelblue") +
  geom_line(aes(y = I_c), color = "black", linetype=2) + 
  geom_line(aes(y = I_u), color = "black") +
  geom_line(aes(y = R_c), color = "darkolivegreen4", linetype=2) +
  geom_line(aes(y = R_u), color = "darkolivegreen4") +
  ggtitle("Process Model")

plot(numeric.sol)

### Export simulated dataset
write_csv(SEIR_sims, "SEIR_sims.csv", append = FALSE)


############# Input Simulated Data into a pomp object ############# 

SEIR_sims %>% 
  select(Week,reports_c,reports_u) %>%
  pomp(
    times="Week",t0=30,
    rprocess=euler(SEIR.compart,delta.t=1/7),
    rinit=SEIR_rinit,
    rmeasure=NULL,
    dmeasure=dmeas,
    accumvars=accumvars,
    statenames=statenames,
    paramnames=paramnames
  ) -> SEIR

############# Weekly Case Data - THE REAL DATA ############# 

# Import dataset for weekly case counts in the wide format

read.csv("weekly_wide.csv", header=TRUE) %>%
  select(Week=MMWRweek, reports_c=PUL, reports_u=WSU) -> WSU_PUL

WSU_PUL %>% as.data.frame() %>% head()

cases.plot <- ggplot(WSU_PUL, aes(x = Week)) + theme_minimal()+ ylab("Case Reports") + xlab("Week")+
  geom_line(aes(y = reports_c, color = 'reports_c')) + 
  geom_line(aes(y = reports_u, color = 'reports_u')) + guides(colour="none") + ggtitle("The Data")

plot(cases.plot)


############# Input Real Data into a pomp object ############# 

WSU_PUL %>%
  pomp(
    times="Week",t0=30,
    rprocess=euler(SEIR.compart,delta.t=1/7),
    rinit=SEIR_rinit,
    rmeasure=rmeas,              # cannot simulate without rmeas
    dmeasure=dmeas,              # only need dmeas for the pfilter
    accumvars=accumvars,
    statenames=statenames,
    paramnames=paramnames
  ) -> SEIR


############# Test the Particle Filter (pf) and compute SE for the pf

foreach (i=1:10, .combine=c) %dopar% {library(pomp)
  SEIR %>% 
    pfilter(Np=5000, params=true.theta)
} -> pf 

logLik(pf) -> ll
logmeanexp(ll,se=TRUE)

plot(pf)


### Check Output

# ess (effective sample size) - number of independent particles and needs to be above 10-100


############# pMCMC Setup ############# 

### Pick variances for the proposal

standard.devs <-c(0.1, 1.0, 0.01, 0.005, 0.1)
params.est <- c("beta_c","beta_u","beta_m","phi","k")
rw.sd <- setNames(standard.devs, params.est)
proposal <- mvn_diag_rw(rw.sd)

# Original Starting Values: standard.devs <-c(0.5, 0.5, 0.001, 0.001, 0.5)

### Set starting values for unknown parameters for n chains

n=5

### Latin Hypercube Sampling

library(lhs)

rows=5
names=c("beta_c", "beta_u", "beta_m", "phi", "k")

X <- randomLHS(rows, length(names))

X[,1] <- qunif(X[,1], 1, 3)
X[,2] <- qunif(X[,2], 20, 30)
X[,3] <- qunif(X[,3], 0.01, 0.2)
X[,4] <- qunif(X[,4], 0.05, 0.15)
X[,5] <- qunif(X[,5], 1, 3)

params <- data.frame(X)
colnames(params) <- names


### Known parameter starting values

alpha=rep(1/3.59,n)
gamma=rep(1/5,n)

N_c=rep(20785,n)
E_c0=rep(0,n)
I_c0=rep(124,n)
R_c0=rep(19,n)

N_u=rep(14254,n)
E_u0=rep(0,n)
I_u0=rep(81,n)
R_u0=rep(15,n)


### Create starting values data frame

theta.start <- data.frame(cbind(params, alpha, gamma, N_c, E_c0, I_c0, R_c0, N_u, E_u0, I_u0, R_u0))
                                

### Run pMCMC

M=101000  # the number of mcmc iterations to run
start_time <- Sys.time()
foreach (theta.start=iter(theta.start,"row"), .inorder=FALSE) %dopar% {library(pomp)
  SEIR %>% 
    pmcmc(
      Nmcmc=M,
      proposal=proposal,
      Np=500,
      params=theta.start
    ) -> pmcmc
  results <- as.data.frame(traces(pmcmc))
} -> results_pmcmc
end_time <- Sys.time()
end_time - start_time    # capture computational time

list_results_pmcmc <- results_pmcmc[c(seq(1:nrow(theta.start)))]
chain.no <- seq(1:n)
iter.no <- seq(1, M+1, by=1)

for(i in seq_along(list_results_pmcmc)){
  list_results_pmcmc[[i]]$chain <- rep(chain.no[i],nrow(list_results_pmcmc[[i]]))
  list_results_pmcmc[[i]]$iter <- iter.no
}

posterior <- do.call(rbind, list_results_pmcmc)

### Export pMCMC Results
write.csv(posterior, "mcmc_results.csv")

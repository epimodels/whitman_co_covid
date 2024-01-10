##### Figure 4 - posterior processing #####

library(tidyverse) # issues with R server
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)
library(coda)
library(reshape2)
library(MCMCglmm)
library(Rmisc)
library(bayestestR) # runs the ci function to get credible intervals
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)
#library(DescTools)

############# Process pMCMC results ############# 

# To read in previously simulated data use this:
posterior <- data.frame(read.csv("mcmc_results.csv", header=TRUE))

post <- posterior %>% 
  dplyr::select(beta_c, beta_u, beta_m, phi, k, chain, iter)

chains <- post %>% 
  dplyr::mutate_at("chain", as.character)

chains.long <- melt(chains, id=c("chain", "iter"))


### Plot the five unprocessed pMCMC chains for quick review

# Quick plot
upc <- ggplot(chains.long, aes(x = iter, y = value, group=chain )) + 
  geom_line(aes(color=chain)) + 
  theme_minimal() +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#999999", "red", "blue")) + 
  facet_wrap(vars(variable), scales='free', ncol = 3)

plot(upc)


# Labelled plot

# Create labels
plot_names <- as_labeller(c('beta_c' = "paste(beta)[C]",
                            'beta_u' = "paste(beta)[U]",
                            'beta_m' = "paste(beta)[M]",
                            'phi' = "paste(phi)",
                            'k' = "k"),label_parsed)

#my_labels <- as_labeller(c(
#"beta_c" = "Beta (C)",
#"beta_u" = "Beta (U)",
#"beta_m" = "Beta (M)",
#"phi"    = "Phi",
#"k"      = "K"))


###### Supplemental S1. Trace Plots Figure ######
upc <- ggplot(data = chains.long,
              aes(x = iter, y = value)) + 
  geom_line(aes(color=chain), alpha = 0.7) + 
  
  guides(color = guide_legend(title = "MCMC Simulation Chains")) +
  
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#999999", "red3", "steelblue")) +
  
  # labels for x-axis
  scale_x_continuous(labels = scales::comma) +
  
  theme_minimal() +
  
  theme(legend.position = c(0.85,0.2),
        legend.text = element_text(hjust = 0.25),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(face = "bold")) + 
  
  ggtitle("Parameter Estimation Random Walk Simulations using pMCMC") +
  xlab("Model Simulations") +  
  ylab("Parameter Values (by scale)") + 
  facet_wrap(vars(variable), scales='free', ncol = 3, labeller = plot_names)

plot(upc)


### Post-process the chains

post.list <- split(post, f = post$chain)  # converts the dataframes back into a list for post-processing
mcmc.list <- mcmc.list(list())

### Fills the list with MCMC objects

for(i in seq_along(post.list)){
  mcmc.list[[i]] <- mcmc(post.list[[i]])
}


### Coda package tools

# Identify burn-in period
raftery.diag(mcmc.list, q=0.025, r=0.005, s=0.95, converge.eps=0.001)

# Convergence Diagnostic - Z-scores should be between 0 and +/-1.96 for convergence
geweke.diag(mcmc.list, frac1=0.1, frac2=0.5)

#geweke.plot(mcmc.list, frac1 = 0.1, frac2 = 0.5, nbins = 20, pvalue = 0.05, auto.layout = TRUE, ask, ...)

# Global convergence tool - gives you output for a point estimate and upper CI. 
# Point estimate needs to be between 1.0 and 1.1 for convergence.
gelman.diag(mcmc.list, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = FALSE)

# Review the autocorrelations for each parameter
c=1 # designate the chain number (1-5)

acf(post.list[[c]]$beta_c, lag.max = 1000)
acf(post.list[[c]]$beta_u, lag.max = 1000) 
acf(post.list[[c]]$beta_m, lag.max = 1000)
acf(post.list[[c]]$phi, lag.max = 1000)
acf(post.list[[c]]$k, lag.max = 1000)

#PlotACF(post.list, lag.max = 10 * log10(length(post.list)), main = NULL, cex = NULL, ...)

M=101000

processed <- window(mcmc.list, start=500, end=M+1, thin=500) 

processed <- data.frame(do.call(rbind, processed))


### Calculate Credible Intervals

ci(processed$beta_c, ci=0.95, method = "HDI") 

ci(processed$beta_u, ci=0.95, method = "HDI")

ci(processed$beta_m, ci=0.95, method = "HDI")

ci(processed$phi, ci=0.95, method = "HDI")

ci(processed$k, ci=0.95, method = "HDI")


### Calculate mode values

processed.long <- melt(processed, id=c("chain", "iter"))

modes <- processed.long %>% 
  dplyr::group_by(variable) %>% 
  dplyr::summarise(mode = posterior.mode(mcmc(value), adjust = 2))

# 0.85 census ascertainment
modes$hdi_low <- c(1.12,19.97,0.00,0.11,0.94)
modes$hdi_high <- c(2.90,39.15,0.15,0.21,2.56)

mode <- as.vector(modes$mode)

modes


# Quick plot for review
P <- ggplot(processed.long, aes(x = value)) + 
  theme_minimal() +
  geom_histogram(aes(y=after_stat(density)), 
                 position="identity", 
                 alpha=0.2, bins=30, colour="black", fill="#999999") +
  geom_density(alpha=.5, fill="#999999", color="black", adjust=2) +
  geom_vline(data=modes, aes(xintercept=mode), color="black", linewidth=1, linetype=2) +
  facet_wrap(vars(variable),labeller = as_labeller(plot_names), scales='free', ncol = 3) + 
  labs(x="Parameter Value",y="Density") +
  theme(panel.spacing=unit(0, "lines"), plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm')) +
  geom_segment(data=modes, aes(x=hdi_low,xend=hdi_high),y=0,yend=0,color="black",linewidth=2,lineend="round") +
  theme(strip.text = element_text(size = 16))

plot(P)


### Final lot with legend for Figure 4

p1 <- ggplot(data = processed.long, aes(x = value)) + 
  
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", bins = 30, 
                 colour = "black", fill = "#999999", alpha = 0.2) +
  
  geom_density(alpha = 0.5, colour = "black", fill = "#999999", adjust = 2) +
  
  geom_vline(data = modes, aes(xintercept = mode), color = "red", linewidth = 1, linetype = 2) + 
  
  geom_segment(data = modes, aes(x=hdi_low, xend=hdi_high), 
               y=0, yend=0, color = "black", linewidth = 1.5, lineend = "round") + 
  
  theme_minimal() +
  
  theme(panel.spacing = unit(0, "lines"), 
        #plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
        strip.text = element_text(size = 16),
        legend.position = c(0.85,0.2),
        legend.text = element_text(hjust = 0.25),
        plot.title = element_text(hjust = 0.5),
        #plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(face = "bold")) + 
  
  ggtitle(label = "Parameter Estimation Posterior Densities using pMCMC",
          subtitle = "") +
  xlab("Parameter Values (by scale)") +  
  ylab("Posterior Density Values (by scale)") + 
  
  facet_wrap(vars(variable), scales = 'free', ncol = 3, labeller = as_labeller(plot_names))

plot(p1)

# Organize parameters for table
Parameter <- c(\beta_C, as.character())

Parameter <- c(expression(beta), expression(phi), parse = TRUE)

Parameter <- c('beta_c' = "paste(beta)[C]",
               'beta_u' = "paste(beta)[U]",
               'beta_m' = "paste(beta)[M]",
               'phi' = "paste(phi)",
               'k' = "k")

mytable <- cbind(Parameter=c('\u03B2c','\u03B2u','\u03B2m','\u03C6','k'), Modes=round(modes$mode, 2))

# Add in the 0.85 census CIs from above
mytable1 <- cbind(mytable, "95% CIs"=c("[1.12, 2.90]","[19.97, 39.15]","[0.00, 0.15]","[0.11, 0.21]","[0.94, 2.56]"))

# Create table grob 
g <- tableGrob(mytable1, rows = NULL, theme = ttheme_minimal(base_size = 7))

g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))

g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))

# Add the table to the figure
ggdraw() +
  draw_plot(p1, 0, 0, 1, 1) +
  draw_plot(g, 0.59, 0.04, 0.5, 0.4)

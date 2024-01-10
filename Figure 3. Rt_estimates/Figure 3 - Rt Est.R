### Rt Estimation of Pullman COVID-19 Epidemics ###

library(dplyr)
library(lubridate)
library(ggplot2)
library(incidence) #Required for EpiEstim package
library(EpiEstim)
library(RColorBrewer)


# Import data
read.csv("WSU_Pullman_wide.csv", header=TRUE) %>%
  select(Report_Date, Week=MMWRweek, Pullman, WSU, Total) -> daily

# Format dates
daily$Report_Date <- mdy(daily$Report_Date)
daily$Report_Date <- as.Date(daily$Report_Date) 

# Filter dates to ensure same length of incidence objects are achieved
daily <- daily %>% filter(Report_Date>='2020-07-26' & Report_Date<='2021-01-02')


### Create incidence objects ###

class(daily) # object should be a dataframe
str(daily)   # check date and integer formats

# Total Cases (reported infections)
incidence_total <- incidence::as.incidence(daily$Total, dates = daily$Report_Date)
plot(incidence_total)

# University Cases (reported infections)
incidence_wsu <- as.incidence(daily$WSU, dates = daily$Report_Date)
plot(incidence_wsu)

# Community Cases (reported infections)
incidence_pullman <- as.incidence(daily$Pullman, dates = daily$Report_Date)
plot(incidence_pullman)


### EpiEstim - Estimate Rt on a sliding weekly window using a parametric serial interval ###

# Now that we have an incidence object we need to supply `EpiEstim` with the SI distribution. 
# As explained above, when we estimate Rt using a parametric SI distribution we only need to 
# supply a mean and standard deviation. 
# We assume that the mean of the SI distribution is 5.19 and that the SD is 1.46 using an estimate 
# found in the literature (Rai B, Shukla A, Dwivedi LK. Estimates of serial interval for COVID-19: 
#                           A systematic review and meta-analysis. Clin Epidemiol Glob Health. 
#                           2021;9:157-161. doi:10.1016/j.cegh.2020.08.007).


# In order to account for the low case count at the start of outbreak, adjusting the start times,
T <- nrow(incidence_total)
t_start <- seq(13, T-6)      # starting at 2 as conditional on the past observations
t_end <- t_start + 6         # adding 6 to get 7-day windows as bounds included in window

R_si_para_Total <- estimate_R(incidence_total, 
                              method = "parametric_si",
                              config = make_config(list(
                                t_start = t_start,
                                t_end = t_end,
                                mean_si=5.19, 
                                std_si=1.46)))

# We supply `estimate_R()` with the incidence object and parameters for the parametric SI and 
# this generates the following 3 panel plot:

plot(R_si_para_Total, legend=FALSE)

# The top panel shows the epidemic curve with the incidence of reported cases over time. 
# The middle panel is the estimated value of R over weekly sliding windows and how it varies 
# over time (Rt), the mean is shown by the solid black line and the 95% credible intervals 
# are represented by the grey shaded area. The final panel is the explored SI distribution 
# that we parameterised and used to estimate Rt.

# You will notice in the "Estimated R" panel that there is a lot of uncertainty around the 
# initial estimates for R~t~ in the outbreak due to very little data. This is why we recommend 
# that you do not try to estimate R~t~ until there have been *at least* 10 cases so that you 
# can be confident that your estimates reflect the data (Cori 2013).

plot(R_si_para_Total, legend=FALSE, "R")


### Estimate Rt for the university and community sub-populations

# University Cases
T <- nrow(incidence_wsu)
w_start <- seq(13, T-6) # starting at 16 as conditional on the past observations relieves the warning
w_end <- w_start + 6 # adding 6 to get 7-day windows as bounds included in window

R_si_para_WSU <- estimate_R(incidence_wsu, 
                            method = "parametric_si",
                            config = make_config(list(
                              t_start = w_start,
                              t_end = w_end,
                              mean_si=5.19,
                              std_si=1.46)))


plot(R_si_para_WSU, legend=FALSE)

plot(R_si_para_WSU, legend=FALSE, "R")


# Community Cases
T <- nrow(incidence_pullman)
p_start <- seq(13, T-6) # starting at 30 as conditional on the past observations relieves the warning
p_end <- p_start + 6 # adding 6 to get 7-day windows as bounds included in window

R_si_para_Pullman <- estimate_R(incidence_pullman, 
                                method = "parametric_si",
                                config = make_config(list(
                                  t_start = p_start,
                                  t_end = p_end,
                                  mean_si=5.19,
                                  std_si=1.46)))

plot(R_si_para_Pullman, legend=FALSE)

plot(R_si_para_Pullman, legend=FALSE, "R")


### -----------------------------------------------------------

# Create a total cases dataset with R estimated and CIs

R_total <- R_si_para_Total$R

R_total <- R_total %>% 
  select(t_start, t_end, mean='Mean(R)', stdev='Std(R)', 
         l_bound='Quantile.0.05(R)', h_bound='Quantile.0.95(R)', 
         median='Median(R)')

######### IMPORTANT NOTE #########  

# Ensure the dataset above contains equal observations as the date values within the incidence object.
# If not, pad with zeros if rows are added.

# Export the dataset for assessment
# write.csv(R_total, "R_total.csv")

##################################

# Import the dataset with correct dimensions (the dataset is included and available for import)
R_total <- read.csv("R_total.csv") %>%
  select(t_start,t_end,mean,stdev,l_bound,h_bound,median)

# Add the dates to the Rt estimates

R_dates <- R_si_para_Total$dates
R_total <- cbind(R_dates, R_total)


### Plot both Rt estimates and bounds for Figure 3a ------------------------------------

ggplot(data = R_total, aes(x = R_dates)) +
  geom_line(aes(y = mean), color = "black", lwd = 1) +
  geom_ribbon(aes(ymin = l_bound, ymax = h_bound), color = "grey", fill = "darkgrey", alpha = 0.6) +
  theme_classic() + 
  theme(legend.position = c(0.9,0.6),
        legend.text = element_text(hjust = 0.25),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(face = "bold")) +
  #axis.text.x = element_text(angle = 45, size = 4.5)) +
  #ggplot2::annotate("text", x = as.numeric(1), y = as.numeric(2), label = "R = 1", color = "red") +
  geom_hline(yintercept=1, linetype = "dashed", lwd = 1.3, color = "darkslategrey") +
  #scale_fill_discrete(name = "Population Type", labels = c("University", "Community")) + 
  ggtitle(label = "Estimation of the Time-varying Reproduction Number (Rt) for SARS-CoV-2 Reported Cases, 
Pullman, WA", 
          subtitle = "Aug-Dec 2020") +
  #ylab(label = "Estimated R values", subtitle = "w/ 95% credible intervals") +
  #xlab(label = "Time Period (MMWR Week 31-53)") +
  labs(x = "Time Period (MMWR Week 31-53)", 
       y = "Estimated R values",
       color = "Population Type")


### -----------------------------------------------------------

# Create university dataset with R estimated and CIs

R_wsu <- R_si_para_WSU$R

R_wsu <- R_wsu %>% 
  select(t_start, t_end, mean='Mean(R)', stdev='Std(R)', 
         l_bound='Quantile.0.05(R)', h_bound='Quantile.0.95(R)', 
         median='Median(R)')

######### IMPORTANT NOTE #########  

# Ensure the dataset above contains equal observations as the date values within the incidence object.
# If not, pad with zeros if rows are added.

# Export the dataset for assessment
# write.csv(R_wsu, "R_wsu.csv")

##################################

# Import the dataset with correct dimensions (the dataset is included and available for import)
R_wsu <- read.csv("R_wsu.csv")


# Create community dataset with R estimated and CIs

R_pul <- R_si_para_Pullman$R

R_pul <- R_pul %>% 
  select(t_start, t_end, mean='Mean(R)', stdev='Std(R)', 
         l_bound='Quantile.0.05(R)', h_bound='Quantile.0.95(R)', 
         median='Median(R)')


######### IMPORTANT NOTE #########  

# Ensure the dataset above contains equal observations as the date values within the incidence object.
# If not, pad with zeros if rows are added.

# Export the dataset for assessment
# write.csv(R_pul, "R_pul.csv")

##################################

# Import the dataset with correct dimensions (the dataset is included and available for import)
R_pul <- read.csv("R_pul.csv")


# Add dates to the Rt estimates

R_dates <- R_si_para_Pullman$dates
R_wsu <- cbind(R_dates, R_wsu)

R_pul <- cbind(R_dates, R_pul)

R_pul %>% 
  select(R_dates=R_dates, mean_p=mean, stdev_p=stdev, 
         l_bound_p=l_bound, h_bound_p=h_bound, median_p=median) -> R_pul


# Combine datasets for plotting
R_data <- left_join(R_wsu, R_pul, by = "R_dates")


### Plot both Rt estimates and bounds for Figure 3b ------------------------------------

ggplot(data = R_data, aes(x = R_dates)) +
  
  geom_line(aes(y = mean, color = "University (w/ 95% CI)"), 
            stat = "identity") +
  geom_ribbon(aes(ymin = l_bound, ymax = h_bound, color = "95% CI (Uni)"), 
              stat = "identity", color = "lightskyblue1", fill = "#377eb8", alpha = 0.5) +
  geom_line(aes(y = mean_p, color = "Community (w/ 95% CI)"), 
            stat = "identity") +
  geom_ribbon(aes(ymin = l_bound_p, ymax = h_bound_p), 
              stat = "identity", color = "lightpink", fill = "#e41a1c", alpha = 0.3) +
   
  theme_classic() + 
  theme(legend.position = c(0.9,0.6),
        legend.text = element_text(hjust = 0.25),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(face = "bold")) +
  
  # add R = 1 line
  geom_hline(yintercept=1, linetype = "dashed", lwd = 1.3, color = "darkslategrey") +

  ggtitle(label = "Estimation of the Time-varying Reproduction Number (Rt) for SARS-CoV-2 Reported Cases 
by Sub-population, Pullman, WA", 
          subtitle = "Aug-Dec 2020") +
  labs(x = "Time Period (MMWR Week 31-53)", 
       y = "Estimated R values",
       color = "Sub-population")

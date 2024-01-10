###### Figure 2. Epidemiological Curves

library(dplyr)
library(lubridate)
library(ggplot2)
library(RColorBrewer)


### Import & Format Weekly Data - Figure 2a

read.csv("weekly_wide.csv", header=TRUE) %>%
  select(Week_date=EpiWeek, Week=MMWRweek, Community=PUL, University=WSU, Total) -> WSU_PUL

WSU_PUL$Week_date <- mdy(WSU_PUL$Week_date)
WSU_PUL$Week_date <- as.Date(WSU_PUL$Week_date)
WSU_PUL$Community <- as.numeric(WSU_PUL$Community)
WSU_PUL$University <- as.numeric(WSU_PUL$University)
WSU_PUL$Total <- as.numeric(WSU_PUL$Total)

### Import & Format Weekly Data - Figure 2b

read.csv("weekly_long_ext.csv", header=TRUE) %>%
  select(Week_date=EpiWeek, Week=MMWRweek, n_cases, Subpopulation=Population) -> long_count_ext

long_count_ext$Week_date <- mdy(long_count_ext$Week_date)
long_count_ext$Week_date <- as.Date(long_count_ext$Week_date)
long_count_ext$n_cases <- as.numeric(long_count_ext$n_cases)
long_count_ext$Subpopulation <- as.factor(long_count_ext$Subpopulation)


### Figure 2a - Epi Curve Total Case Population

epi_tot <- ggplot(data = WSU_PUL, aes(x = Week_date)) + # x-axis is epiweek (as class Date)
  
  geom_bar(aes(y = Total),          # y-axis height in the weekly case counts
           stat = "identity", 
           width = 4, 
           fill = "grey35", 
           alpha = 0.6) +  
  
  # labels for x-axis
  scale_x_date(
    date_breaks = "1 month",        # labels every 2 months 
    date_minor_breaks = "1 month",  # gridlines every month
    date_labels = '%b\n%Y') +       # labeled by month with year below
  
  # Adjust color palette if preferred (uses RColorBrewer package)
  #scale_fill_brewer(palette = "Set1") + 
  
  theme_classic() +
  
  theme(
    legend.position = c(0.85,0.7),
    legend.text = element_text(hjust = 0.25),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line = element_line(color = 'black'),
    axis.title = element_text(face = "bold")) +    
  
  labs(
    title    = "Total Weekly Reported SARS-CoV-2 Cases, Pullman, WA",
    subtitle = "Aug-Dec 2020",
    x        = "MMWR Week",
    y        = "Weekly Cases Reported",
    fill     = "Total Cases")

plot(epi_tot)


### Figure 2b - Epi curve by sub-population

epi_sub <- ggplot(data = long_count_ext, aes(x = Week_date,)) + # x-axis is epiweek (as class Date)
  
  geom_histogram(
    mapping = aes(
      y = n_cases,                # y-axis height in the weekly case counts
      group = Subpopulation,     
      fill = Subpopulation),
    stat = "identity",
    position = "dodge",
    alpha = 0.6) +     
  
  # labels for x-axis
  scale_x_date(
    date_breaks = "1 month",         # labels every 2 months 
    date_minor_breaks = "1 month",   # gridlines every month
    date_labels = '%b\n%Y') +        # labeled by month with year below
  
  # Choose color palette (uses RColorBrewer package)
  scale_fill_brewer(palette = "Set1") + 
  
  theme_classic() +
  
  theme(
    legend.position = c(0.85,0.7),
    legend.text = element_text(hjust = 0.25),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line = element_line(color = 'black'),
    axis.title = element_text(face = "bold")) + 
  
  labs(
    title    = "Weekly Reported SARS-CoV-2 Cases by Sub-population, Pullman, WA",
    subtitle = "Aug-Dec 2020",
    x        = "MMWR Weeks 31-53",
    y        = "Weekly Cases Reported",
    fill     = "Sub-population")

plot(epi_sub)


# devtools::install_github("HighlanderLab/RandPedigreePCA", subdir = "randPedPCA",
#                         ref="v0.9.3", build_vignettes = T)

rm(list = ls())
library(randPedPCA)
library(readr)
library(tidyverse)

Labrador <- read_csv(file = "/Users/roscraddock/Documents/University\ of\ Edinburgh/HighlanderLab - Data_2025/LabPed2025.csv")
Labrador<- Labrador[order(Labrador$gen),]

# Add coat colour? Four official categories: Black, Chocolate, Liver, Yellow then unknown.
Labrador$Coat_Colour <- Labrador$Colour
labcolours <- c("BLACK", "YELLOW", "LIVER", "CHOCOLATE")
Labrador <- Labrador %>% 
  transform(Coat_Colour = ifelse(Coat_Colour %in% labcolours, Coat_Colour, "Unknown/nonStandard"))
Labrador <- Labrador %>%
  transform(Coat_Colour = ifelse(Coat_Colour == "BLACK", "Black", Coat_Colour))
Labrador <- Labrador %>%
  transform(Coat_Colour = ifelse(Coat_Colour == "YELLOW", "Yellow", Coat_Colour))
Labrador <- Labrador %>%
  transform(Coat_Colour = ifelse(Coat_Colour == "LIVER", "Chocolate", Coat_Colour))
Labrador <- Labrador %>%
  transform(Coat_Colour = ifelse(Coat_Colour == "CHOCOLATE", "Chocolate", Coat_Colour))

# Create Pedigree Object
ped <- pedigree(Labrador$fatherID,
                Labrador$motherID,
                Labrador$id)

# Obtain centred estimate from L inverse using Hutch++
li <- sparse2spam(getLInv(ped)) # generate L inverse and convert to spam format
set.seed(123345) # set random seed as Hutch++ uses random numbers
hp <- hutchpp(li,num_queries = 100, center = T) # estimate

# Run randPedPCA
pc <- rppca(ped, center=T, totVar = hp)
summary(pc)

# Collect proportion of variances explained by PC1 and PC2
var <- pc[["varProps"]]*100
# Labels for x and y axis
pc1 <- paste("PC1 (",round(var[1], 2), "%)", sep = "")
pc2 <- paste("PC2 (",round(var[2], 2), "%)", sep = "")

# Collect PC in t
t <- pc[["x"]]
# Plot PC1 and PC2, labelled by coat colour
ggplot(data = t, aes(x = t[,1], y = t[,2], colour = Labrador$Coat_Colour)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Labrador Pedigree",
       x = pc1,
       y = pc2,
       colour = "Coat Colour") +
  scale_colour_manual(values = c("Black" = "black",
                                 "Yellow" = "yellow",
                                 "Chocolate" = "brown",
                                 "Unknown/nonStandard" = "grey"),
                      breaks = c("Black", "Yellow", "Chocolate", "Unknown/nonStandard"),
                      labels = c("Black", "Yellow", "Chocolate", "Non Standard/ Unknown")) +
                        theme(
                          axis.title.x = element_text(size = 22),  # Increase x-axis label text size
                          axis.title.y = element_text(size = 22),  # Increase y-axis label text size
                          plot.title = element_text(size = 24),    # Increase plot title text size
                          legend.title = element_text(size = 22),  # Increase legend title text size
                          legend.text = element_text(size = 20),
                          axis.text.x = element_text(size = 20),   # Increase x-axis tick labels text size
                          axis.text.y = element_text(size = 20))    # Increase y-axis tick labels text size

# Animate PCA
library(gifski)
library(gganimate)

Labrador$year <- as.integer(Labrador$year)

# Collect frames of the PCA plot labelled by coat colour for each year of birth. Total: 70 frames
gapplot <- ggplot(data = t, aes(x = t[,1], y = t[,2], colour = Labrador$Coat_Colour)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Labrador Pedigree, Year: {frame_time}",
       x = pc1,
       y = pc2,
       colour = "Coat Colour") +
  scale_colour_manual(values = c("Black" = "black",
                                 "Yellow" = "yellow",
                                 "Chocolate" = "brown",
                                 "Unknown/nonStandard" = "grey"),
                      breaks = c("Black", "Yellow", "Chocolate", "Unknown/nonStandard"),
                      labels = c("Black", "Yellow", "Chocolate", "Non Standard/ Unknown")) +
  theme(
    axis.title.x = element_text(size = 22),  # Increase x-axis label text size
    axis.title.y = element_text(size = 22),  # Increase y-axis label text size
    plot.title = element_text(size = 24),    # Increase plot title text size
    legend.title = element_text(size = 22),  # Increase legend title text size
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 20),   # Increase x-axis tick labels text size
    axis.text.y = element_text(size = 20))  +  # Increase y-axis tick labels text size
  transition_time(Labrador$year) +
  ease_aes('linear')

# Save the animated plot as a gif
animate(gapplot, height = 600, width = 800, nframes = 70, fps = 5, renderer=gifski_renderer("LabradorPCA2.gif"))


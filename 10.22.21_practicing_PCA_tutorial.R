#practicing PCA
#source: https://datavizpyr.com/how-to-make-pca-plot-with-r/

library(tidyverse)
library(broom)
library(palmerpenguins)

penguins <- penguins %>%
  drop_na() %>%
  select(-year)

#this dataset has several variables: species, island, bill length, bill depth, flipper length, body mass, sex
#(and year, which we just removed)

#"We use select() to select numerical variables in penguins’s data, apply scale() and then do PCA with prcomp() function."

pca_fit <- penguins %>%
  select(where(is.numeric)) %>%
  scale() %>%
  prcomp()

#Here is a quick summary of the PCA. We can see that the first two principal components explain 88% of the variation in the data.
summary(pca_fit)

#"Our PCA results do not contain any “meta” information and the original data. 
#We will use broom’s augment() function to add the original data to pca results."
pca_fit %>%
  augment(penguins)

#"rename the PCs columns to PC1, PC2,…and so on."
pca_fit %>%
  augment(penguins) %>%
  rename_at(vars(starts_with(".fitted")),
            list(~str_replace(.,".fitted","")))

#"Now we have the data ready for making a PCA plot, in this example a scatter plot between the first two Principal Components. 
#Since we have the original data handy, we can color the data points by species variable and change the shape by sex variable."
pca_fit %>%
  augment(penguins) %>%
  rename_at(vars(starts_with(".fitted")),
            list(~str_replace(.,".fitted",""))) %>%
  ggplot(aes(x=PC1, 
             y=PC2,
             color=species,
             shape=sex))+
  geom_point()

#okay, so what does this mean for my own data?


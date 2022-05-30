#rm(list = ls())
#set.seed(123)

library(mcr)
library(ggplot2)
library(MASS)
library(tidyverse)
library(ggpmisc)
library(GGally)
library(MASS)
library(caret)
library(skimr)
library(patchwork)




##--------------------------------------- Set Standard Deviation ---------------------------------------
sims <- c(0.1,0.2,0.3)


## Simulating Dataset ---------------------------------------------------------------------------------

sigma <- rbind(c(1,0,0), c(0,1,0), c(0,0,1))  ## No Corr, Unit Variance
mu <- c(3, 5, 1)  # create the mean vector


#OUTER LOOP START -------------------------------------------------------------------------------------

nested_var_list = list()

for (k in 69:100) {  #try with 20, set.seed


V_df <- as.data.frame(mvrnorm(n=300, mu=mu, Sigma=sigma)) # Generate multivariate normal dist.


#Extra variables
V4 <- rpois(n = 300, lambda = 4)         #poisson
V5 <- rnorm(n = 300, mean = 25, sd = 4)  #normal


#LOOP START --------------------------------------------------------------------------------------------

  datalist_var = list()

  for (j in 1:length(sims)) {


    Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
    e <- rnorm(n=300, mean = 0, sd=sims[j]*Y) # Error generated as a function of response
    Y <- Y + e

    df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe

    ones <- rep(1,dim(df)[1])  # Generating column of ones matching intercept for suff stat calculations
    X <- cbind(ones,df)
    X <- as.matrix(X[,-2])

    fref <- lm(Y ~ . , data = df) # Reference Model
    mean(fref$resid^2)


    ## MCR - Binary Search --------------------------------------------------------------------------------

    p <- dim(X)[2]
    MCRminus <- MCRplus <- c(1:p)


    for (i in 1:p){

      suff_stat <- get_suff_stats_lm(Y, X, p1 = i)
      objects <- precompute_mcr_objects_and_functions(X, Y, p1 = i, model_class_loss = 'linear_mse')


      fref_loss <- get_e0_lm(model = fref$coefficients, suff_stats = suff_stat)

      MCR2 <- get_empirical_MCR(eps = fref_loss + 0.05*fref_loss,objects) #.05 arbitrary
      MCRminus[i] <- MCR2$range[1]
      MCRplus[i] <- MCR2$range[2]
    }

    mcr_table <- data.frame(id = c('V1','V2','V3','V4','V5'),MCRminus=round(MCRminus[-1],4),MCRplus=round(MCRplus[-1],4))
    datalist_var[[j]] <- mcr_table

  }

nested_var_list[[k]] <- datalist_var

}

# Save an object to a file
saveRDS(nested_var_list, file = "sim_heteroscedasticty_output1.rds")
# Restore the object
nested_var_list1 <-readRDS(file = "sim_heteroscedasticty_output1.rds")


#Data Manipulation
#-----------------------------------------------------------------------------
data_var_1 <- mapply(data.frame, lapply(nested_var_list1,'[[',1) ,SIMPLIFY = F)
data_var_2 <- mapply(data.frame, lapply(nested_var_list1,'[[',2) ,SIMPLIFY = F)
data_var_3 <- mapply(data.frame, lapply(nested_var_list1,'[[',3) ,SIMPLIFY = F)
#-----------------------------------------------------------------------------


wide1 <- reshape::merge_all(data_var_1)
long1 <- wide1 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p1 <- ggplot(long1, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,9) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Heteroscedasticity Level 1') + xlab('Variables') + ylab('Model Class Reliance') +
  scale_fill_brewer(palette="Dark2") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                               plot.title = element_text(hjust = 0.5))

wide2 <- reshape::merge_all(data_var_2)
long2 <- wide2 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p2 <- ggplot(long2, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,9) + geom_boxplot(width=0.95)  + geom_hline(yintercept = 1) +
  ggtitle('Heteroscedasticity Level 2') + xlab('Variables') + ylab(NULL) +
  scale_fill_brewer(palette="Dark2") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                               plot.title = element_text(hjust = 0.5))

wide3 <- reshape::merge_all(data_var_3)
long3 <- wide3 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p3 <- ggplot(long3, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,9) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Heteroscedasticity Level 3') + xlab('Variables') + ylab(NULL) +
  scale_fill_brewer(palette="Dark2") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                               plot.title = element_text(hjust = 0.5))

combined <- p1 + p2 + p3  & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")
ggsave('figures/heteroscedasticity.png', combined)


#Heteroscedasticity Plots(Funnels-----------------------

sigma <- rbind(c(1,0,0), c(0,1,0), c(0,0,1))  ## No Corr, Unit Variance
mu <- c(3, 5, 1)  # create the mean vector
V_df <- as.data.frame(mvrnorm(n=300, mu=mu, Sigma=sigma)) # Generate multivariate normal dist.

V4 <- rpois(n = 300, lambda = 4)         #poisson
V5 <- rnorm(n = 300, mean = 25, sd = 4)  #normal


#1
Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
e <- rnorm(n=300, mean = 0, sd=0.1*Y) # Error generated as a function of response
Y <- Y + e
df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe

fref1 <- lm(Y ~ . , data = df) # Reference Model
mean(fref1$resid^2)

h1 <- ggplot(mapping = aes(y = resid(fref1),
                     x = fitted(fref1))) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point() +
  theme_minimal() + ylim(-12,12)+
  xlab('Fitted Values') + ylab('Residuals') +ggtitle('Heteroscedasticity Level 1')
#2
Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
e <- rnorm(n=300, mean = 0, sd=0.2*Y) # Error generated as a function of response
Y <- Y + e
df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe

fref2 <- lm(Y ~ . , data = df) # Reference Model
mean(fref2$resid^2)


h2 <- ggplot(mapping = aes(y = resid(fref2),
                           x = fitted(fref2))) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point() +
  theme_minimal() + ylim(-12,12)  +
  xlab('Fitted Values') + ylab('Residuals') +ggtitle('Heteroscedasticity Level 2')

#3
Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
e <- rnorm(n=300, mean = 0, sd=0.3*Y) # Error generated as a function of response
Y <- Y + e
df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe

fref3 <- lm(Y ~ . , data = df) # Reference Model
mean(fref3$resid^2)


h3 <- ggplot(mapping = aes(y = resid(fref3),
                           x = fitted(fref3))) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point() +
  theme_minimal() + ylim(-12,12) +
  xlab('Fitted Values') + ylab('Residuals') +ggtitle('Heteroscedasticity Level 3')

combined <- h1 + h2 + h3
combined + plot_layout(ncol = 3)


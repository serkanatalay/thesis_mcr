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


#------------------------------ Set Mean and Covariances below ---------------------------------------

mu<-c(2, 5, 1)  # create the mean vector

sigma <- rbind(c(1,0,0), c(0,1, 0), c(0,0,1))   # No Collinearity

sample_n = c(100,300,500)

title_list <- list('n=100','n=300','n=500')


#OUTER LOOP START -------------------------------------------------------------------------------------
nested_n_list = list()

for (k in 1:100) {  #100 repetitions

#LOOP START --------------------------------------------------------------------------------------------

  datalist_n = list()

  for (j in 1:(length(sample_n))){

    V_df<-as.data.frame(mvrnorm(n=sample_n[j], mu=mu, Sigma=sigma))

    #Extra variables
    V4 <- rpois(n = sample_n[j], lambda = 4)         #poisson
    V5 <- rnorm(n = sample_n[j], mean = 25, sd = 4)  #normal

    Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
    e <- rnorm(n=sample_n[j], mean = 0, sd=1.5)
    Y <- Y + e

    df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe
    head(df,n = 3)

    ones <- rep(1,dim(df)[1])  # Generating column of ones matching intercept for suff stat calculations
    X <- cbind(ones,df)
    X <- as.matrix(X[,-2])

    fref <- lm(Y ~ . , data = df) # Reference Model


    p <- dim(X)[2]
    MCRminus <- MCRplus <- c(1:p)

    for (i in 1:p){

      suff_stat <- get_suff_stats_lm(Y, X, p1 = i)
      objects <- precompute_mcr_objects_and_functions(X, Y, p1 = i, model_class_loss = 'linear_mse')
      #}

      fref_loss <- get_e0_lm(model = fref$coefficients, suff_stats = suff_stat)

      MCR2 <- get_empirical_MCR(eps = fref_loss + .05*fref_loss,objects) #.25 arbitrary
      MCRminus[i] <- MCR2$range[1]
      MCRplus[i] <- MCR2$range[2]
    }

    mcr_table <- data.frame(id = c('V1','V2','V3','V4','V5'),MCRminus=round(MCRminus[-1],4),MCRplus=round(MCRplus[-1],4))
    datalist_n[[j]] <- mcr_table

  }

  nested_n_list[[k]] <- datalist_n

}


# Save to a file
saveRDS(nested_n_list, file = "sim_sample_size_output.rds")
# Restore the object
nested_n_list <-readRDS(file = "sim_sample_size_output.rds")

#Data Manipulation
#-----------------------------------------------------------------------------
data_n_1 <- mapply(data.frame, lapply(nested_n_list,'[[',1) ,SIMPLIFY = F)
data_n_2 <- mapply(data.frame, lapply(nested_n_list,'[[',2) ,SIMPLIFY = F)
data_n_3 <- mapply(data.frame, lapply(nested_n_list,'[[',3) ,SIMPLIFY = F)
#-----------------------------------------------------------------------------



#Graphs
##-----------------------------------------------------------------------------
wide1 <- reshape::merge_all(data_n_1)
long1 <- wide1 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p1 <- ggplot(long1, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,7.90) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Sample Size = 100') + xlab('Variables') + ylab('Model Class Reliance') +
  scale_fill_brewer(palette="Set3") + theme_minimal() + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

wide2 <- reshape::merge_all(data_n_2)
long2 <- wide2 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p2 <- ggplot(long2, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,7.90) + geom_boxplot(width=0.95)  + geom_hline(yintercept = 1) +
  ggtitle('Sample Size = 300') + xlab('Variables') +  ylab(NULL) +
  scale_fill_brewer(palette="Set3") + theme_minimal() + theme(legend.position="none",panel.border=element_rect(fill=NA),plot.title = element_text(hjust = 0.5))

wide3 <- reshape::merge_all(data_n_3)
long3 <- wide3 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p3 <- ggplot(long3, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,7.90) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Sample Size = 500') + xlab('Variables') + ylab(NULL) +
  scale_fill_brewer(palette="Set3") + theme_minimal() + theme(legend.position="none",panel.border=element_rect(fill=NA),plot.title = element_text(hjust = 0.5))

combined <- p1 + p2 + p3 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")
ggsave('figures/sample_size.png', combined)

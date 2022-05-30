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

aList <- list(rbind(c(1, 0  ,0), c( 0,  1, 0), c(0,0,1)),         # level0 No-Collinearity
              rbind(c(1,-0.2,0), c(-0.2,1, 0), c(0,0,1)),         # level1 Collinearity
              rbind(c(1,-0.4,0), c(-0.4,1, 0), c(0,0,1)),         # level2 Collinearity
              rbind(c(1,-0.6,0), c(-0.6,1, 0), c(0,0,1)),         # level3 Collinearity
              rbind(c(1,-0.8,0), c(-0.8,1, 0), c(0,0,1)))         # level4 Collinearity

title_list <- list('level0 No-Collinearity','level1 Collinearity','level2 Collinearity','level3 Collinearity','level4 Collinearity')

#----------------------------------------------------------------------------------------------------

nested_list = list()

for (k in 1:100) {  #100 repetitions


  mu<-c(2, 5, 1)  # create the mean vector

  #Extra variables
  V4 <- rpois(n = 300, lambda = 4)         #poisson
  V5 <- rnorm(n = 300, mean = 25, sd = 4)  #normal


  #LOOP START --------------------------------------------------------------------------------------------

  datalist_col = list()

  for (j in 1:length(aList)){

    V_df<-as.data.frame(mvrnorm(n=300, mu=mu, Sigma=aList[[j]]))
    Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
    e <- rnorm(n=300, mean = 0, sd=1.5)
    Y <- Y + e

    df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe

    ones <- rep(1,dim(df)[1])  # Generating column of ones matching intercept for suff stat calculations
    X <- cbind(ones,df)
    X <- as.matrix(X[,-2])

    fref <- lm(Y ~ . , data = df) # Reference Model


    p <- dim(X)[2]
    MCRminus <- MCRplus <- c(1:p)

    for (i in 1:(p+1)){

      # last row is permutation of 1 and 2
      if (i == (p+1)){
        suff_stat <- get_suff_stats_lm(Y, X, p1 = c(2,3))
        objects <- precompute_mcr_objects_and_functions(X, Y, p1 = c(2,3), model_class_loss = 'linear_mse')
      }

      # permute each variable individually
      else{
      suff_stat <- get_suff_stats_lm(Y, X, p1 = i)
      objects <- precompute_mcr_objects_and_functions(X, Y, p1 = i, model_class_loss = 'linear_mse')
      }

      fref_loss <- get_e0_lm(model = fref$coefficients, suff_stats = suff_stat)

      MCR2 <- get_empirical_MCR(eps = fref_loss + .05*fref_loss,objects) #.5 arbitrary
      MCRminus[i] <- MCR2$range[1]
      MCRplus[i] <- MCR2$range[2]
    }

    summary(fref)
    mcr_table <- data.frame(id = c('V1','V2','V3','V4','V5','V1 & V2'),MCRminus=round(MCRminus[-1],4),MCRplus=round(MCRplus[-1],4))

    datalist_col[[j]] <- mcr_table


   }

  nested_list[[k]] <- datalist_col

}

# Save to a file
saveRDS(nested_list, file = "sim_coll_levels_output.rds")
# Restore the object
nested_list <-readRDS(file = "sim_coll_levels_output.rds")



#Data Manipulation
#-----------------------------------------------------------------------------
data_col_1 <- mapply(data.frame, lapply(nested_list,'[[',1) ,SIMPLIFY = F)
data_col_2 <- mapply(data.frame, lapply(nested_list,'[[',2) ,SIMPLIFY = F)
data_col_3 <- mapply(data.frame, lapply(nested_list,'[[',3) ,SIMPLIFY = F)
data_col_4 <- mapply(data.frame, lapply(nested_list,'[[',4) ,SIMPLIFY = F)
data_col_5 <- mapply(data.frame, lapply(nested_list,'[[',5) ,SIMPLIFY = F)




#Graphs
##-----------------------------------------------------------------------------
wide1 <- reshape::merge_all(data_col_1)
long1 <- wide1 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p1 <- ggplot(long1, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,8.5) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('No - Multicollinearity') + xlab('Variables') + ylab('Model Class Reliance') +
  scale_fill_brewer(palette="Pastel1") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

wide2 <- reshape::merge_all(data_col_2)
long2 <- wide2 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p2 <- ggplot(long2, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,8.5) + geom_boxplot(width=0.95)  + geom_hline(yintercept = 1) +
  ggtitle('Multicollinearity - 0.2') + xlab('Variables') + ylab('Model Class Reliance') +
  scale_fill_brewer(palette="Pastel1") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),plot.title = element_text(hjust = 0.5))

wide3 <- reshape::merge_all(data_col_3)
long3 <- wide3 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p3 <- ggplot(long3, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,8.5) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Multicollinearity - 0.4') + xlab('Variables') + ylab(NULL) +
  scale_fill_brewer(palette="Pastel1") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

wide4 <- reshape::merge_all(data_col_4)
long4 <- wide4 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p4 <- ggplot(long4, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,8.5) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Multicollinearity - 0.6') + xlab('Variables') + ylab('Model Class Reliance') +
  scale_fill_brewer(palette="Pastel1") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

wide5 <- reshape::merge_all(data_col_5)
long5 <- wide5 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p5 <- ggplot(long5, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,8.5) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Multicollinearity - 0.8') + xlab('Variables') + ylab(NULL) +
  scale_fill_brewer(palette="Pastel1") + theme_minimal(base_size=11.5) + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

combined <- p1 + p2 + p3 +  p4 + p5 & theme(legend.position = "bottom")

layout <- c(
  area(1, 2, 1, 3),
  area(2, 1, 2, 2),
  area(2, 3, 2, 4),
  area(3, 1, 3, 2),
  area(3, 3, 3, 4))
combined <- combined + plot_layout(guides = "collect",design = layout)
ggsave('figures/multicol.png', combined)

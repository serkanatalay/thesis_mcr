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



#------------------------------ Rashomon Error Thresholds below ---------------------------------------

threshold_lev <- c(0.01,0.025,0.05)


#OUTER LOOP START -------------------------------------------------------------------------------------
nested_tre_list = list()

for (k in 1:100) {  #100 repetitions

#------------------------------ Generating the dataset ------------------------------------------------

  mu <- c(2, 5, 1)  # create the mean vector

  sigma <- rbind(c(1,0,0), c(0 ,1, 0), c(0,0,1))  ## NO COLLINEARITY

  V_df<-as.data.frame(mvrnorm(n=300, mu=mu, Sigma=sigma))

  #Extra variables
  V4 <- rpois(n = 300, lambda = 4)         #poisson
  V5 <- rnorm(n = 300, mean = 25, sd = 4)  #normal

  Y <- 2*V_df$V1 + 1*V_df$V3 + 1*V4 # Target variable without error
  e <- rnorm(n=300, mean = 0, sd=1.5)
  Y <- Y + e

  df<- data.frame(cbind(Y,V_df,V4,V5)) # Creating Dataframe
  head(df,n = 3)

  ones <- rep(1,dim(df)[1])  # Generating column of ones matching intercept for suff stat calculations
  X <- cbind(ones,df)
  X <- as.matrix(X[,-2])


  #LOOP START --------------------------------------------------------------------------------------------

  datalist_tre = list()

  for (j in 1:(length(threshold_lev))){


    fref <- lm(Y ~ . , data = df) # Reference Model

    p <- dim(X)[2]
    MCRminus <- MCRplus <- c(1:p)

    for (i in 1:p){

      suff_stat <- get_suff_stats_lm(Y, X, p1 = i)
      objects <- precompute_mcr_objects_and_functions(X, Y, p1 = i, model_class_loss = 'linear_mse')
      #}

      fref_loss <- get_e0_lm(model = fref$coefficients, suff_stats = suff_stat)

      MCR2 <- get_empirical_MCR(eps = fref_loss + threshold_lev[j]*fref_loss,objects)
      MCRminus[i] <- MCR2$range[1]
      MCRplus[i] <- MCR2$range[2]
    }

    mcr_table <- data.frame(id = c('V1','V2','V3','V4','V5'),MCRminus=round(MCRminus[-1],4),MCRplus=round(MCRplus[-1],4))
    datalist_tre[[j]] <- mcr_table

  }

nested_tre_list[[k]] <- datalist_tre

}

# Save to a file
saveRDS(nested_tre_list, file = "sim_err_threshold_output.rds")
# Restore the object
nested_tre_list <-readRDS(file = "sim_err_threshold_output.rds")


##-----------------------------------------------------------------------------
data_tre_1 <- mapply(data.frame, lapply(nested_tre_list,'[[',1) ,SIMPLIFY = F)
data_tre_2 <- mapply(data.frame, lapply(nested_tre_list,'[[',2) ,SIMPLIFY = F)
data_tre_3 <- mapply(data.frame, lapply(nested_tre_list,'[[',3) ,SIMPLIFY = F)

#
# sd01 <- data_tre_1 %>%  dplyr::bind_rows() %>% group_by(id) %>%
#   mutate(MCRrange = MCRplus-MCRminus) %>%  summarise_all(list(mean = mean, sd = sd)) %>% #mean, sd
#   select(MCRrange_sd)
# sd02 <- data_tre_2 %>%  dplyr::bind_rows() %>% group_by(id) %>%
#   mutate(MCRrange = MCRplus-MCRminus) %>% summarise_all(list(mean = mean, sd = sd)) %>% #mean, sd
#   select(MCRrange_sd)
# sd03 <- data_tre_3 %>%  dplyr::bind_rows() %>% group_by(id) %>%
#   mutate(MCRrange = MCRplus-MCRminus) %>% summarise_all(list(mean = mean, sd = sd)) %>% #mean, sd
#   select(MCRrange_sd)
#
# xtable(bind_cols(sd01,sd02,sd03))

##-----------------------------------------------------------------------------
#MCR Ranges t.test
#threshold1
V1_tr01 <- data_tre_1 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>%  filter(id=='V1') %>% select(MCRrange)
V2_tr01 <- data_tre_1 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V2') %>% select(MCRrange)
V3_tr01 <- data_tre_1 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V3') %>% select(MCRrange)
V4_tr01 <- data_tre_1 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V4') %>% select(MCRrange)
V5_tr01 <- data_tre_1 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V5') %>% select(MCRrange)
#threshold2
V1_tr025 <- data_tre_2 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>%  filter(id=='V1') %>% select(MCRrange)
V2_tr025 <- data_tre_2 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V2') %>% select(MCRrange)
V3_tr025 <- data_tre_2 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V3') %>% select(MCRrange)
V4_tr025 <- data_tre_2 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V4') %>% select(MCRrange)
V5_tr025 <- data_tre_2 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V5') %>% select(MCRrange)
#threshold3
V1_tr05 <- data_tre_3 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>%  filter(id=='V1') %>% select(MCRrange)
V2_tr05 <- data_tre_3 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V2') %>% select(MCRrange)
V3_tr05 <- data_tre_3 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V3') %>% select(MCRrange)
V4_tr05 <- data_tre_3 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V4') %>% select(MCRrange)
V5_tr05 <- data_tre_3 %>%  dplyr::bind_rows() %>%
  mutate(MCRrange = MCRplus/MCRminus) %>% filter(id=='V5') %>% select(MCRrange)

t.test(V1_tr01, V1_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V2_tr01, V2_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V3_tr01, V3_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V4_tr01, V4_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V5_tr01, V5_tr05, alternative = "two.sided", var.equal = FALSE)
##-----------------------------------------------------------------------------
t.test(V1_tr025, V1_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V2_tr025, V2_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V3_tr025, V3_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V4_tr025, V4_tr05, alternative = "two.sided", var.equal = FALSE)
t.test(V5_tr025, V5_tr05, alternative = "two.sided", var.equal = FALSE)
##-----------------------------------------------------------------------------
t.test(V1_tr025, V1_tr01, alternative = "two.sided", var.equal = FALSE)
t.test(V2_tr025, V2_tr01, alternative = "two.sided", var.equal = FALSE)
t.test(V3_tr025, V3_tr01, alternative = "two.sided", var.equal = FALSE)
t.test(V4_tr025, V4_tr01, alternative = "two.sided", var.equal = FALSE)
t.test(V5_tr025, V5_tr01, alternative = "two.sided", var.equal = FALSE)
##-----------------------------------------------------------------------------


#Graphs
##-----------------------------------------------------------------------------


wide1 <- reshape::merge_all(data_tre_1)
long1 <- wide1 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p1 <- ggplot(long1, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,6.5) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Rashomon Threshold - 0.01') + xlab('Variables') + ylab('Model Class Reliance') +
  scale_fill_brewer(palette="Spectral") + theme_minimal() + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

wide2 <- reshape::merge_all(data_tre_2)
long2 <- wide2 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p2 <- ggplot(long2, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,6.5) + geom_boxplot(width=0.95)  + geom_hline(yintercept = 1) +
  ggtitle('Rashomon Threshold - 0.025') + xlab('Variables') +  ylab(NULL) +
  scale_fill_brewer(palette="Spectral") + theme_minimal() + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

wide3 <- reshape::merge_all(data_tre_3)
long3 <- wide3 %>% pivot_longer(cols=c(MCRminus,MCRplus),names_to = 'type',values_to = 'MCR')
p3 <- ggplot(long3, aes(x=id, y=MCR, fill=type )) +  ylim(0.5,6.5) + geom_boxplot(width=0.95) + geom_hline(yintercept = 1) +
  ggtitle('Rashomon Threshold - 0.05') + xlab('Variables') + ylab(NULL) +
  scale_fill_brewer(palette="Spectral") + theme_minimal() + theme(legend.position="none",panel.border=element_rect(fill=NA),
                                                                  plot.title = element_text(hjust = 0.5))

combined <- p1 + p2 + p3 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

combined <- combined + plot_layout(guides = "collect", design = layout)
ggsave('figures/threshold.png', combined)

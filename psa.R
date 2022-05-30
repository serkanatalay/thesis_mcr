library(tidyr)
library(mcr)
library(tidyverse)
library(ggpmisc)
library(GGally)
library(skimr)
library(corrplot)
library(caret)
library(MASS)
library(xlsx)



## Obtain MCR ranges for all of the variables

dat1 <- read.table("prostate.data", header = TRUE)
dat1 <- dat1[,-10]

ggpair_plot_psa <- ggpairs(dat1[,c(9,1,5,2,3,4)])
ggsave('figures/ggpair', ggpair_plot_psa)

ones <- rep(1,dim(dat1)[1])
y <- dat1$lpsa
X <- matrix(as.numeric(unlist(dat1[,c(-7,-9)])),nrow=nrow(dat1[,c(-7,-9)]))
X <- as.matrix(cbind(ones, X))


fref <- lm(lpsa ~ . - gleason , data = dat1) # fit model without intercept; Run this to obtain the coefficients
summary(fref)


## Bootstrap starts ------------------------------------------------------------

boots <- 1:dim(dat1)[1]
boot_list <- list()

for (b in 1:100){

  e <- sample(boots, size = length(boots), replace = TRUE)
  y.boot <- y[e]
  X.boot <- X[e,]

  ## MCR - Binary Search ------------------------------------------------------------

  p <- dim(X)[2]
  MCRminus <- MCRplus <- c(1:p)

  for (i in 1:p){

    # permute each variable individually

      suff_stat <- get_suff_stats_lm(y.boot, X.boot, p1 = i)
      objects <- precompute_mcr_objects_and_functions(X.boot, y.boot, p1 = i, model_class_loss = 'linear_mse')


    fref_loss <- get_e0_lm(model = fref$coefficients, suff_stats = suff_stat)

    MCR2 <- get_empirical_MCR(eps = fref_loss + .025*fref_loss,objects) #.0.25 arbitrary
    MCRminus[i] <- MCR2$range[1]
    MCRplus[i] <- MCR2$range[2]
  }
  mcr_table <- data.frame(id = c('ones',colnames(dat1[,c(-7,-9)])),MCRminus=round(MCRminus,4),MCRplus=round(MCRplus,4))
  boot_list[[b]] <- mcr_table
}


# Save an object to a file
saveRDS(boot_list, file = "boot_list_psa.rds")
# Restore the object
boot_list_psa <-readRDS(file = "boot_list_psa.rds")


#Quantiles
qdf <- reshape::merge_all(boot_list_psa) %>%
  group_by(id) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = 'id', values_from=c('MCRminus','MCRplus')) %>%
  select(-contains(c("row","ones"))) %>%
  summarise_all(list(q025 =~quantile(.,probs=0.025,na.rm=TRUE), q975=~quantile(.,probs=0.975,na.rm=TRUE)))


#Datalist to Dataframe
data_var_1 <- mapply(data.frame, lapply(boot_list_psa,'[[',1) ,SIMPLIFY = F)
data_var_2 <- mapply(data.frame, lapply(boot_list_psa,'[[',2) ,SIMPLIFY = F)
data_var_3 <- mapply(data.frame, lapply(boot_list_psa,'[[',3) ,SIMPLIFY = F)


summary_auto <- reshape::merge_all(boot_list) %>% group_by(id) %>% summarize_all(list(mean=mean,sd=sd,median=median))
summary_auto

## Plot ---------------------------------------------------------------------------
ggplot(mcr_table, aes(x=id))+
  geom_linerange(aes(ymin=MCRminus,ymax=MCRplus),linetype=2,color="gray45")+
  geom_point(aes(y=MCRminus,color="MCR-"),size=3.25)+
  geom_point(aes(y=MCRplus,color="MCR+"),size=3.25)+
  scale_y_continuous(breaks = seq(0, 15, by = 1)) +
  geom_hline(yintercept=1,linetype="dashed") +
  theme_bw() + labs(y = "MCR Range", x = 'Variables') +
  scale_colour_manual(values=c("red","dodgerblue2")) +
  geom_table_npc(data = mcr_table, label = list(mcr_table[order(mcr_table$id),])
                 ,npcx = 0.02, npcy = 0.95, hjust = 0, vjust = 1)

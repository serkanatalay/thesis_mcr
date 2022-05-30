library(tidyverse)
library(ggpmisc)
library(GGally)
library(skimr)
library(corrplot)
library(mcr)
library(caret)
library(MASS)
library(xlsx)
library(caret)




## preprocessing -------------------------------------------------------------------

auto <- read.table('autompg.data',header = FALSE)
auto <- auto[c(-9)] #getting rid of car model names

auto$V4 <- as.numeric(auto$V4) #convert to numeric
auto$V8 <- as.factor(auto$V8)

auto <-drop_na(auto)

dummy <- dummyVars(" ~ .", data=auto)
auto <- data.frame(predict(dummy, newdata = auto))

column_names <- c('mpg','cylinders','displacement','horsepower','weight','acceleration','model_year','origin_us','origin_eu','origin_jp')
colnames(auto) <- column_names





## Obtain MCR ranges for all of the variables -------------------------------------
fref <- lm(mpg ~ .-1 , data = auto)
summary(fref)

ggpair_plot_auto <- ggpairs(auto[c(1,3,5,7,8)]) #pairplot significant variables
#ggsave('figures/ggpair_auto', ggpair_plot_auto)



y <- auto$mpg
#ones <- rep(1,dim(auto)[1])
#X <- auto[-1] #add column of ones for suff_stat calculations (intercept)
X <- matrix(as.numeric(unlist(auto[-1])),nrow=nrow(auto[-1]))
#X <- as.matrix(X)
head(X)
head(auto[-1])



## Bootstrap starts ------------------------------------------------------------

boots <- 1:dim(auto)[1]
boot_list <- list()

for (b in 1:100){

  e <- sample(boots, size = length(boots), replace = TRUE)
  y.boot <- y[e]
  X.boot <- X[e,]

## MCR - Binary Search ------------------------------------------------------------

  p <- dim(X)[2]
  MCRminus <- MCRplus <- c(1:(p-2))

  for (i in 1:(p-2)){

    # last row is permutation of 1 and 2
    if (i == (p-2)){
      suff_stat <- get_suff_stats_lm(y.boot, X.boot, p1 = c(7,8,9))
      objects <- precompute_mcr_objects_and_functions(X.boot, y.boot, p1 = c(7,8,9), model_class_loss = 'linear_mse')}

    # permute each variable individually
    else{
      suff_stat <- get_suff_stats_lm(y.boot, X.boot, p1 = i)
      objects <- precompute_mcr_objects_and_functions(X.boot, y.boot, p1 = i, model_class_loss = 'linear_mse')
      }

      fref_loss <- get_e0_lm(model = fref$coefficients, suff_stats = suff_stat)

      MCR2 <- get_empirical_MCR(eps = fref_loss + .05*fref_loss,objects) #.0.05 arbitrary
      MCRminus[i] <- MCR2$range[1]
      MCRplus[i] <- MCR2$range[2]
    }

  mcr_table <- data.frame(id = c(column_names[-c(1,8,9,10)],'origin'),MCRminus=round(MCRminus,4),MCRplus=round(MCRplus,4))
  boot_list[[b]] <- mcr_table

  }

# Save an object to a file
saveRDS(boot_list, file = "boot_list.rds")
# Restore the object
boot_list <-readRDS(file = "boot_list.rds")


#Quantiles
qdfauto <- reshape::merge_all(boot_list) %>%
  group_by(id) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = 'id', values_from=c('MCRminus','MCRplus')) %>%
  select(everything(), -contains(c("row","ones"))) %>%
  summarise_all(list(q025 =~quantile(.,probs=0.025,na.rm=TRUE), q975=~quantile(.,probs=0.975,na.rm=TRUE)))

#write.xlsx(as.data.frame(qdfauto),file='qdfauto.xlsx','Sheet2',row.names = F)


#Datalist to Dataframe
data_var_1 <- mapply(data.frame, lapply(boot_list,'[[',1) ,SIMPLIFY = F)
data_var_2 <- mapply(data.frame, lapply(boot_list,'[[',2) ,SIMPLIFY = F)
data_var_3 <- mapply(data.frame, lapply(boot_list,'[[',3) ,SIMPLIFY = F)


summary_auto <- reshape::merge_all(boot_list) %>% group_by(id) %>% summarize_all(list(mean=mean,sd=sd,median=median))
summary_auto

#quantile(.95)

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
    ,npcx = 0.05, npcy = 0.95, hjust = 0, vjust = 1)

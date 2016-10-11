---
  title: "RT_variability"
---
  
  install.packages("tidyr"); install.packages("plyr"); install.packages("dplyr"); install.packages("ggplot2")
install.packages("latticeExtra"); install.packages("stats"); install.packages("zoo"); install.packages("Rmisc");
install.packages("tibble"); install.packages("stringr"); install.packages("data.table");


suppressMessages(library(tidyr));suppressMessages(library(plyr));suppressMessages(library(dplyr)); 
suppressMessages(library(ggplot2)); suppressMessages(library(latticeExtra)); suppressMessages(library(stats));
suppressMessages(library(zoo)); suppressMessages(library(Rmisc)); suppressMessages(library(tibble)); 
suppressMessages(library(stringr)); suppressMessages(library(data.table));

### Load beh data 

#study beh
basedir = '/Volumes/emmt/'
subs = read.csv(paste0(basedir,'sublists/subs_eeg_beh_subset_n25.txt'), header=F)
subs = subs$V1
f_list <- unlist(lapply(subs, 
                        function(s) {
                          s = str_pad(s, 3, pad='0')
                          paste0(basedir,'data/beh/s',s,'/eMMT',s,'_study_concat_mem.csv')
                        }))

length(f_list)

raw_rt <- do.call('rbind',
                  lapply(f_list,
                         function(x) {
                           sub_d <- fread(x, header = T, sep=',') %>% as.data.frame()
                           sub_d <- rownames_to_column(sub_d, var = "Trial")
                         }
                  )
) 

#test beh 
f_list_test <- unlist(lapply(subs, 
                             function(s) {
                               s = str_pad(s, 3, pad='0')
                               paste0(basedir,'data/beh/s',s,'/eMMT',s,'_test_concat_mem.csv')
                             }))

length(f_list_test)

raw_rt_test <- do.call('rbind',
                       lapply(f_list_test,
                              function(x) {
                                sub_d <- fread(x, header = T, sep=',') %>% as.data.frame()
                                sub_d <- rownames_to_column(sub_d, var = "Trial")
                              }
                       )
) %>%
  filter(mem != "noresponse")

d_rt <- raw_rt
d_rt_test <- raw_rt_test

#calculate controlled study RT (subtracting trial RT by average RT for that word) 
d_rt <- d_rt %>% 
  group_by(word) %>%
  mutate(mean_RT = mean(RT)) %>% 
  mutate(controlled_study_RT = RT - mean_RT) %>%
  ungroup() %>%
  group_by(subjNum) %>%
  mutate(controlled_study_RT_zscore = abs(as.numeric(scale(controlled_study_RT, center = TRUE, 
                                                           scale = TRUE))))

#calculate mem_accuracy for each participant 
#load data
word_acc = read.csv('/Volumes/emmt/analysis/mturk_accuracy/word_char_acc.csv')
#calculate accuracy rate
d_rt <- word_acc %>% 
  mutate(category = ifelse(p_human > 0.75, "human",
                           ifelse(p_human < 0.25, "nonhuman", "either"))) %>%
  select(word, category) %>% 
  right_join(d_rt, by = "word", copy=FALSE) %>% 
  select(subjNum, word, category, resp, controlled_study_RT_zscore, onset, mem) %>%
  mutate(accuracy = ifelse(resp == "j" & category == "human", 1,
                           ifelse(resp == "j" & category == "nonhuman", 0,
                                  ifelse(resp == "k" & category == "human", 0,
                                         ifelse(resp == "k" & category == "nonhuman", 1, 1))))) %>%
  mutate(interpolate_zscore = ifelse(accuracy == 0, NA , controlled_study_RT_zscore)) %>%
  group_by(subjNum) %>%
  mutate(interpolate_zscore = ((na.locf(interpolate_zscore, na.rm = FALSE)) + 
                                 (na.locf(interpolate_zscore, na.rm = FALSE, fromLast = TRUE)) / 2))

#calculate smoothed zscore using controlled RT
d_rt <- d_rt %>% 
  group_by(subjNum) %>% 
  mutate(smoothed_onset = ksmooth(onset, interpolate_zscore, "normal", bandwidth = 16)[[1]],
         smoothed_zscore = ksmooth(onset, interpolate_zscore, "normal", bandwidth = 16)[[2]])

#calculate zone using smoothed_RT_var
d_rt <- d_rt %>%
  group_by(subjNum) %>%
  mutate(controlledRT.zone = ifelse(smoothed_zscore> median(smoothed_zscore), "out",
                                    ifelse(smoothed_zscore < median(smoothed_zscore), "in", "")))

#calculate accuracy data 
accuracy.data <- d_rt %>%
  group_by(subjNum, controlledRT.zone) %>%
  mutate(accuracy_rate = sum(accuracy)/n()) %>%
  select(subjNum, accuracy_rate, controlledRT.zone) %>%
  unique(incomparables = FALSE) %>%
  filter(controlledRT.zone != "NA" | controlledRT.zone != "")

#plot accuracy and controlledRT
accuracy.data$subjNum <- as.factor(accuracy.data$subjNum)
ggplot(accuracy.data, aes(x=controlledRT.zone, y=accuracy_rate)) + 
  geom_line(aes(colour=subjNum, group=subjNum)) +
  geom_point(aes(colour=subjNum), size=3)

ggplot(accuracy.data, aes(x=controlledRT.zone, y=accuracy_rate)) +
  geom_jitter(position=position_jitter(0.2), cex=1.2) + 
  stat_summary(fun.data=data_summary, color="blue")

### calculate d' data (controlled for memorability)

#only select words from IQR in terms of memorability
word_hits_rt <- d_rt %>% 
  group_by(word, mem) %>% 
  summarise(item_memorability=n()) %>% 
  filter(mem == 'hi_hit') %>% 
  select(word, item_memorability) %>%
  mutate(mem_bin = ntile(item_memorability, 4)) %>%
  filter(mem_bin == 2 || mem_bin == 3)

mem_words <- d_rt %>%
  filter(word %in% word_hits_rt$word) %>%
  select(subjNum, word, controlledRT.zone)

#based on word in IQR, calculate dprime for in the zone, and out of the zone words
d_rt_test <- d_rt_test %>% 
  left_join(mem_words, by = c("word", "subjNum"), copy = FALSE)

in_zone <- d_rt_test %>%
  filter(controlledRT.zone == "in" | is.na(controlledRT.zone)) %>%
  group_by(subjNum) %>%
  mutate(dprime = individual.dprime(mem))

out_zone <- d_rt_test %>%
  filter(controlledRT.zone == "out" | is.na(controlledRT.zone)) %>%
  group_by(subjNum) %>%
  mutate(dprime = individual.dprime(mem))

zone_dprime <- in_zone %>%
  bind_rows(out_zone) %>%
  select(subjNum, controlledRT.zone, dprime) %>%
  filter(controlledRT.zone == "in" | controlledRT.zone == "out") %>%
  unique()

##plot zone and mem (dprime)
zone_dprime$subjNum <- as.factor(zone_dprime$subjNum)
ggplot(zone_dprime, aes(x=controlledRT.zone, y=dprime)) + 
  geom_line(aes(colour=subjNum, group=subjNum)) +
  geom_point(aes(colour=subjNum), size=3)

ggplot(zone_dprime, aes(x=controlledRT.zone, y=dprime)) +
  geom_jitter(position=position_jitter(0.2), cex=1.2) + 
  stat_summary(fun.data=data_summary, color="blue")


##FUNCTIONS  
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#calculate quantile
Quantile <- function (v, p) {
  v = sort(v)
  m = 0
  n = length(v)
  j = floor((n * p) + m)
  g = (n * p) + m - j
  y = ifelse (g == 0, 0, 1)
  ((1 - y) * v[j]) + (y * v[j+1])
}

#calculate individual dprime
individual.dprime <- function(currentsub){
  hit <- sum(currentsub == "hi_hit")
  miss <- sum(currentsub == "hi_miss")
  fa <- sum(currentsub == "hi_fa")
  cr <- sum(currentsub == "hi_cr")
  hit <- (hit/(hit+miss))
  fa <- (fa/(fa+cr))
  return (qnorm(hit) - qnorm(fa))
}


#plotting individual subject smoothed zscore
ggplot(filter(d, subjNum == 12), aes(x=smoothed_onset, y=smoothed_zscore)) + geom_line()

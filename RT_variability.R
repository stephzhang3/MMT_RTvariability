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

basedir = '/Volumes/emmt/'
subs = read.csv(paste0(basedir,'sublists/subs_beh.txt'), header=F)
subs = subs$V1
f_list <- unlist(lapply(subs, 
                        function(s) {
                          s = str_pad(s, 3, pad='0')
                          paste0(basedir,'data/beh/s',s,'/eMMT',s,'_study_concat_mem.csv')
                        }))

length(f_list)

d <- do.call('rbind',
             lapply(f_list,
                    function(x) {
                      sub_d <- fread(x, header = T, sep=',') %>% as.data.frame()
                      sub_d <- rownames_to_column(sub_d, var = "Trial")
                    }
             )
)

#calculate controlled study RT (subtracting trial RT by average RT for that word) 
d <- d %>% 
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
d <- word_acc %>% 
  mutate(category = ifelse(p_human > 0.75, "human",
                           ifelse(p_human < 0.25, "nonhuman", "either"))) %>%
  select(word, category) %>% 
  right_join(zone.sub.data, by = "word", copy=FALSE) %>% 
  select(subjNum, word, category, resp, controlled_study_RT_zscore, onset) %>%
  mutate(accuracy = ifelse(resp == "j" & category == "human", 1,
                           ifelse(resp == "j" & category == "nonhuman", 0,
                                  ifelse(resp == "k" & category == "human", 0,
                                         ifelse(resp == "k" & category == "nonhuman", 1, 1))))) %>%
  mutate(interpolate_zscore = ifelse(accuracy == 0, NA , controlled_study_RT_zscore)) %>%
  group_by(subjNum) %>%
  mutate(interpolate_zscore = ((na.locf(interpolate_zscore, na.rm = FALSE)) + 
                                 (na.locf(interpolate_zscore, na.rm = FALSE, fromLast = TRUE)) / 2))

#calculate smoothed zscore using controlled RT
d <- d %>% 
  group_by(subjNum) %>% 
  mutate(smoothed_onset = ksmooth(onset, interpolate_zscore, "normal", bandwidth = 16)[[1]],
         smoothed_zscore = ksmooth(onset, interpolate_zscore, "normal", bandwidth = 16)[[2]])

#calculate zone using smoothed_RT_var
d <- d %>%
  group_by(subjNum) %>%
  mutate(controlledRT.zone = ifelse(smoothed_zscore> median(smoothed_zscore), "out",
                                    ifelse(smoothed_zscore < median(smoothed_zscore), "in", "")))

accuracy.data <- d %>%
  group_by(subjNum, controlledRT.zone) %>%
  mutate(accuracy_rate = sum(accuracy)/n()) %>%
  select(subjNum, accuracy_rate, controlledRT.zone) %>%
  unique(incomparables = FALSE) %>%
  filter(controlledRT.zone != "NA") %>%
  filter(controlledRT.zone != "")

#plot accuracy and controlledRT
accuracy.data$subjNum <- as.factor(d$subjNum)
ggplot(accuracy.data, aes(x=controlledRT.zone, y=accuracy_rate)) + 
  geom_line(aes(colour=subjNum, group=subjNum)) +
  geom_point(aes(colour=subjNum), size=3)

ggplot(accuracy.data, aes(x=controlledRT.zone, y=accuracy_rate)) +
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


#plotting individual subject smoothed zscore
ggplot(filter(d, subjNum == 12), aes(x=smoothed_onset, y=smoothed_zscore)) + geom_line()

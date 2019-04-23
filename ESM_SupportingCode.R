#'---
#' title: "Supporting code for comment on Jones et al."
#' author: "Larkin, Buck, Fieberg, Galatowitsch"
#' output: 
#'    html_document:
#'       toc: true
#'       theme: default
#'       toc_depth: 2
#'       toc_float:
#'           collapsed: false
#'---

#' ### Load libraries and data
#+ message=FALSE, warning=FALSE
options(width=160)
library(tidyr)
library(dplyr)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringr)
library(geepack)
library(metafor)
library(lme4)
library(car)

#' Data and function from JEA 'Code for ProcB paper.R'
load("ProcB paper.RData") 
#' - "latdata" = full dataset from JEA Dryad (5,142 records)
#' - "outIQR" = function to remove outliers based on IQR (lines 44-51 in JEA code)
#' 
#' # 1. Issues with Recovery Completeness 
#' 
#' #### ("Response Ratio" in JEA's dataset)
#' 
#' ## 1.1 Sign reversal
#' 
#' Confirming that sign reversal based on relationship between End and Goal, not Start 
#' and Goal. 
#' 
#' Re-creating values reported by JEA in 'latdata' by adding 0.05 when Goal or End = 0 
#' and reversing sign for cases where End > Goal
#+ fig.height = 4
latdata2 <- 
  latdata %>%
  mutate(adjEnd = ifelse(End==0 | Goal==0, End + 0.05, End), 
         adjStart = ifelse(End==0 | Goal==0, Start + 0.05, Start),
         adjGoal = ifelse(End==0 | Goal==0, Goal + 0.05, Goal), 
         calcResponseRatio = ifelse(adjEnd < adjGoal, log(adjEnd/adjGoal), -1*log(adjEnd/adjGoal)))
par(mar = c(5,5,1,1))
plot(ResponseRatio ~ calcResponseRatio, data = latdata2)  

#' Subsetting data w/ complete Start, Goal, and End data to identify restoration outcomes 
#' found in the data (since progress cannot be evaluated without knowing starting condition)  
latdata_sge <-
  latdata %>%
  filter(complete.cases(Start, Goal, End)) %>%
  filter(Start != Goal) # eliminating cases where Start = Goal (restoration moot)
save(latdata_sge, file = "latdata_sge.RData")

#' Classifying restoration outcomes based on Start-End-Goal relationship
latdata_sge <-
  latdata_sge %>%
  mutate(goal_direction = ifelse(Goal < Start, "decrease", "increase")) %>%
  mutate(outcome = ifelse(goal_direction == "increase", 
                          ifelse(End < Start, "worse", 
                                 ifelse(End == Start, "nochange",
                                        ifelse(End < Goal, "incomplete",
                                               ifelse(End > Goal, "exceeded", "match")))),
                          ifelse(goal_direction == "decrease", 
                                 ifelse(End > Start, "worse", 
                                        ifelse(End == Start, "nochange",
                                               ifelse(End > Goal, "incomplete",
                                                      ifelse(End < Goal, "exceeded", "match")))), "error")))

#' ## 1.2 False inferences
#' 
#' Perfect matching of the goal was always correctly classified by JEA, but there were both false 
#' exceeds (few) and false incompletes (many), i.e., cases where restoration should have been 
#' interpreted as exceeding the goal but was instead classified as incomplete recovery 
#' due to issues with sign reversal. 
#' 
#' Secondarily, cases where restoration worsened the condition, did not change it, 
#' or yielded incomplete recovery could not be distinguished.
outcomes <-
  latdata_sge %>%
  group_by(goal_direction, outcome) %>%
  summarise(n=n(), n_zero = length(ResponseRatio[ResponseRatio == 0]),
            n_neg = length(ResponseRatio[ResponseRatio < 0]), 
            n_pos = length(ResponseRatio[ResponseRatio > 0])) %>%
  mutate(RR_should_be = ifelse(outcome %in% c("worse", "no_change", "incomplete"), "neg", 
                               ifelse(outcome == "match", "0", "pos")))
outcomes <- outcomes[,c(1:3,7,4:6)]
outcomes

outcomes.prop <-
  outcomes[,c(2:3)] %>%
  group_by(outcome) %>%
  summarise(n = sum(n), n.all = sum(outcomes$n), prop = n/n.all)
outcomes.prop

#' When restoration should have been classified as incomplete recovery, what percentage
#' of these cases were classified as exceeding goal?
outcomes.summ <-
  outcomes[,-c(1)] %>%
  group_by(outcome, RR_should_be) %>%
  summarise(n = sum(n), sum_neg = sum(n_neg), sum_zero = sum(n_zero), sum_pos = sum(n_pos))
outcomes.summ
false.pos <- (outcomes.summ[outcomes.summ$outcome=="incomplete",'sum_pos']/
                outcomes.summ[outcomes.summ$outcome=="incomplete",'n'])*100; false.pos

#' When restoration should have been classified as exceeding goal, what percentage
#' of these cases were classified as incomplete recovery?
false.neg <- (outcomes.summ[outcomes.summ$outcome=="exceeded",'sum_neg']/
                outcomes.summ[outcomes.summ$outcome=="exceeded",'n'])*100; false.neg

#' What cases did JEA identify as exceeding goal (RR > 0) and were these classified correctly? 
#' 
#' RR > 0 when negative values for both End and Goal. 'outcome' is how these should have been 
#' classified based on Start condition
exceed.goal <-
  latdata_sge %>%
  filter(ResponseRatio > 0) %>%
  group_by(goal_direction, outcome) %>%
  summarise(n=n(), min.Goal = min(Goal), max.Goal = max(Goal), 
            min.End = min(End), max.End = max(End))  
exceed.goal

#' Of the cases JEA identified as exceeding goal, <1/3 were accurately classified
exceed.goal.sum <-
  exceed.goal %>%
  group_by(outcome) %>%
  summarise(n = sum(n))
exceed.goal.sum$n[exceed.goal.sum$outcome=="exceeded"]/sum(exceed.goal.sum$n)  

#' ## 1.3 Summary figures
#' 
#' Identifying example cases w/ same Ends and Goals, different outcomes (see R file for code)
#+ echo=FALSE
matches <- find.matches(latdata_sge[, c("Goal", "End")], latdata_sge[, c("Goal", "End")], 
                        tol = c(0.0, 0.0), maxmatch = 100)
matches <- matches$matches
matches <- matches[rowSums(matches) > matches[,1],]
matches <- matches[!duplicated(matches),]
rownames(matches) <- matches[,1]
matches <- melt(matches[,-1])
matches <- matches[!matches$value==0, -2]
names(matches) <- c("Row", "MatchingRow")
matches$Goal <- latdata_sge$Goal[matches$Row]
matches$End <- latdata_sge$End[matches$Row]
matches$Start1 <- latdata_sge$Start[matches$Row]
matches$Start2 <- latdata_sge$Start[matches$MatchingRow]
matches$goal_direction1 <- latdata_sge$goal_direction[matches$Row]
matches$goal_direction2 <- latdata_sge$goal_direction[matches$MatchingRow]
matches$outcome1 <- latdata_sge$outcome[matches$Row]
matches$outcome2 <- latdata_sge$outcome[matches$MatchingRow]
matches$ResponseRatio1 <- latdata_sge$ResponseRatio[matches$Row]
matches$ResponseRatio2 <- latdata_sge$ResponseRatio[matches$MatchingRow]
matches <- matches[matches$Start1 != matches$Start2,]
matches <- matches[matches$outcome1 != "match",]
matches <- matches[matches$outcome1 != matches$outcome2,]

#' Example cases scored as exceeding goals (see R file for code)
#+ echo=FALSE
exceed.goal.exs <-
  latdata_sge %>%
  filter(ResponseRatio > 0) %>%
  select(ID, Citation, Goal, End, Start, goal_direction, outcome, ResponseRatio)
exceed.goal.exs <- latdata_sge[latdata_sge$ID %in% c(14, 3254, 4975),]

#' Choosing subset of these for plotting (see R file for code)
#+ echo=FALSE
exs.dat <- rbind(exceed.goal.exs, latdata_sge[c(1163, 1164, 1256, 590, 147, 163, 2066, 2540),])
exs.dat$ID <- as.character(exs.dat$ID)
exs.dat$row.ID <- row.names(exs.dat)
exs.dat$row.ID <- factor(exs.dat$row.ID, levels = exs.dat$row.ID)
exs.dat$Citation <- as.character(exs.dat$Citation)
exs.dat1 <- exs.dat[-c(1:3),]
exs.dat2 <- exs.dat[c(1:3),]
# Fixing formatting of citation names
exs.dat2[1,]$Citation <- "Riedinger-Whitmore et al. 2005"
exs.dat2[3,]$Citation <- "Waddington & Warner 2001"

#' Plot illustrating mismatches between actual restoration outcomes and 
#' JEA recovery completeness (see R file for code)
#+ fig.height = 4
#+ echo=FALSE
margins <- unit(c(0.25,0.25,0.25,0.25), "cm")
width <- 18
p1a <- ggplot(exs.dat1, aes(row.ID)) +
  geom_point(aes(x = Start, y = row.ID), size = 4, shape = 21, fill = "gray60") +
  geom_segment(aes(x = Start, xend = End, y = row.ID, yend = row.ID), arrow = arrow(angle = 30, length = unit(0.1, "inches"), type = "closed")) +
  geom_point(aes(x = Goal, y = row.ID), size = 4, shape = 21, color = "red", fill = "white") +
  geom_point(aes(x = Goal, y = row.ID), size = 2.5, shape = 4, color = "red") +
  xlab("") +
  ylab("Citation") + 
  scale_y_discrete(labels=str_wrap(exs.dat1$Citation, width = width)) +
  coord_cartesian(ylim = c(1, 8.75), xlim = c(-0.145, 0.535)) +
  geom_text(aes(x = Start+0.01, y = row.ID), label = "Start", data = exs.dat1[c(8),], 
            hjust = 0, vjust = -1.5, size = 3, fontface = 2, angle = 40) +
  geom_text(aes(x = End+0.01, y = row.ID), label = "End", data = exs.dat1[c(8),], 
            hjust = 0, vjust = -1.5, size = 3, fontface = 2, angle = 40) +
  geom_text(aes(x = Goal+0.01, y = row.ID), label = "Goal", data = exs.dat1[c(8),], 
            hjust = 0, vjust = -1.5, size = 3, fontface = 2, angle = 40) +
  geom_text(aes(x = c(exs.dat1$End[1:3], exs.dat1$Start[4], exs.dat1$End[5], exs.dat1$Start[6], 
                      exs.dat1$Goal[7:8]), y = row.ID, 
                label = paste(c("Worsened", "Exceeded goal", "Worsened", "Incomplete", 
                                "Exceeded goal", "No change", "No change", "Incomplete"), 
                              " ", c(rep("(",7),"(Recovery Completeness = "), 
                              sprintf("%0.2f", round(ResponseRatio, digits = 2)), ")", sep = "")), 
            data = exs.dat1, hjust = c(1.05, 1.05, -0.1, -0.15, 1.05, 1.15, -0.15, -0.06), 
            vjust = 0.5, size = 3, fontface = 1) +
  theme(plot.margin = unit(c(0.25,0.25,0,0.25), "cm"))
p1a

#' Plot showing examples where Recovery Completeness > 1 but goal not 
#' exceeded (see R file for code)
#+ fig.height = 2
#+ echo=FALSE
p1b <- ggplot(exs.dat2, aes(row.ID)) +
  geom_point(aes(x = Start, y = row.ID), size = 4, shape = 21, fill = "gray60") +
  geom_segment(aes(x = Start, xend = End, y = row.ID, yend = row.ID), arrow = arrow(angle = 30, length = unit(0.1, "inches"), type = "closed")) +
  geom_point(aes(x = Goal, y = row.ID), size = 4, shape = 21, color = "red", fill = "white") +
  geom_point(aes(x = Goal, y = row.ID), size = 2.5, shape = 4, color = "red") +
  xlab("Ecological measurement") +
  ylab("") + 
  scale_y_discrete(labels=str_wrap(exs.dat2$Citation, width = width)) +
  coord_cartesian(xlim = c(-24, 0)) +
  geom_text(aes(x = c(exs.dat2$End[1], exs.dat2$Start[2], exs.dat2$Goal[3]), y = row.ID,
                label = paste(c("Worsened", "Incomplete", "Incomplete"), " ", c(rep("(",3)), 
                              sprintf("%0.2f", round(ResponseRatio, digits = 2)), ")", sep = "")),
            data = exs.dat2, hjust = c(-0.1, 1.15, 1.15), vjust = 0.4, size = 3, fontface = 1) +
  theme(plot.margin = unit(c(0,0.25,0.25,0.25), "cm"))
p1b

#' Plot illustrating aggregate effects of sign reversal issue on reported
#' restoration outcomes
#' 
#' Outcomes classified based on Start-End-Goal relationship (see R file for code)
#+ echo=FALSE
plot.dat <-
  outcomes[,-1] %>%
  group_by(outcome, RR_should_be) %>%
  summarise(n = sum(n), n_zero = sum(n_zero), n_neg = sum(n_neg), n_pos = sum(n_pos))
plot.dat <- plot.dat[,c(1,3)]
plot.dat$Baseline <- "Start condition" 

#' Outcomes classified based on End-Goal relationship (see R file for code)
#+ echo=FALSE
end.dat <- 
  outcomes[,-c(1:4)] %>%
  gather(key = 'outcome', value = 'n', n_zero:n_pos) %>%  
  group_by(outcome) %>%
  summarise(n = sum(n), Baseline = "End condition\n(JEA method)")
end.dat$outcome <- c("incomplete", "exceeded", "match")

plot.dat <- bind_rows(plot.dat, end.dat)
plot.dat$outcome <- factor(plot.dat$outcome, levels(factor(plot.dat$outcome))[c(5,4,2,3,1)])
plot.dat$Baseline <- factor(plot.dat$Baseline, levels(factor(plot.dat$Baseline))[c(2,1)])

#+ fig.height = 4
#+ echo=FALSE
p2 <- ggplot(plot.dat, aes(x = outcome, y = n, fill = Baseline)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  scale_x_discrete(labels=c("Worsened", "No change", "Incomplete\nrecovery", 
                            "Matched\ngoal", "Exceeded\ngoal")) +
  xlab("Restoration outcomes") +
  ylab("Number of cases") + 
  scale_fill_manual(values = c("gray60", "black")) +
  labs(fill = "Baseline used") + theme(legend.position = c(0.15, 0.7), 
                                       plot.margin = margins)
p2

#' # 2. Issues with Recovery Rate
#' 
#' #### ("Resilience" in JEA's dataset)
#' 
#' ## 2.1 Complete resilience data 
#' 
#' Removing records without resilience values
lat.rs <- latdata[!is.na(latdata$asinhResilience),]
lat.rs <- lat.rs[is.finite(lat.rs$asinhResilience),]

#' Summarizing full resilience data
lat.rs.summ <- 
  lat.rs %>%
  group_by(PA) %>%
  summarise(n = n(), mean_rs = mean(asinhResilience), se_rs = sd(asinhResilience)/sqrt(n))
lat.rs.summ

#' ## 2.2 Data following intended outlier removal
#' 
#' Thresholds for 'asinhResilience' identified using IQR method described in JEA text 
#' differ from cutoffs used in JEA code (line 105), which were -3.304 and 6.608. 
#' Here JEA's 'outIQR' function has been modified to return IQR thresholds. 
outIQR2 <- function(vec) {
  quants = quantile(vec)
  newIQR = 1.5*IQR(vec)
  upper = quants[4] + newIQR
  lower = quants[2] - newIQR
  return(c(lower, upper))
}
outIQR2(lat.rs$asinhResilience)
outlier.lo <- as.numeric(outIQR2(lat.rs$asinhResilience)[1])
outlier.hi <- as.numeric(outIQR2(lat.rs$asinhResilience)[2])

#' Resilience data with outliers removed as intended
#' 
lat.rs.or <- lat.rs[outIQR(lat.rs$asinhResilience),]
lat.rs.or.summ <-   
  lat.rs.or %>%
  filter(asinhResilience > outlier.lo & asinhResilience < outlier.hi) %>%
  group_by(PA) %>%
  summarise(n = n(), mean_rs = mean(asinhResilience), se_rs = sd(asinhResilience)/sqrt(n))
lat.rs.or.summ

#' ## 2.3 Data following implemented outlier removal
#' 
#' Thresholds used were more liberal than IQR method, i.e., discarded more data
lat.rs.or.manual <-   
  lat.rs %>%
  filter(asinhResilience < 6.608 & asinhResilience > -3.304)

lat.rs.or.manual.summ <-   
  lat.rs.or.manual %>%
  group_by(PA) %>%
  summarise(n = n(), mean_rs = mean(asinhResilience), se_rs = sd(asinhResilience)/sqrt(n))
lat.rs.or.manual.summ

#' Size of dataset for full data
sum(lat.rs.summ$n)
#' Size of dataset with outlier removal by IQR
sum(lat.rs.or.summ$n)
#' Size of dataset with errant outlier removal
sum(lat.rs.or.manual.summ$n)
#' Percent reduction in dataset (# of cases) from full data to errant outlier removal
(sum(lat.rs.summ$n) - sum(lat.rs.or.manual.summ$n))/sum(lat.rs.summ$n) * 100
#' Reduction in numbers of cases and citations from intended to errant outlier removal
sum(lat.rs.or.summ$n) - sum(lat.rs.or.manual.summ$n)
length(unique(lat.rs.or$Citation)) - length(unique(lat.rs.or.manual$Citation))
#' Were the discarded data cases where the authors had measured recovery over hours or days?
# summary of full data
summary(lat.rs$TimeSince)
# summary of discarded data
lat.discarded <- lat.rs[is.na(match(lat.rs$ID, lat.rs.or.manual$ID)),]
summary(lat.discarded$TimeSince)
# Proportion w/ TimeSince<1 year
length(lat.discarded$TimeSince[lat.discarded$TimeSince < 1])/length(lat.discarded$TimeSince)
# Minimum discarded time since (in days)
summary(lat.discarded$TimeSince)[1]*365

#' ## 2.4 Summary figures

#' Density plot showing distribution of resilience data for active vs. passive restoration
#' and outlier removal thresholds as JEA intended (dashed lines) and as implemented 
#' (dotted lines/darker grey shading) (see R file for code)
#' 
#+ fig.height = 4
#+ echo = FALSE
p3 <- ggplot(lat.rs, aes(x=asinhResilience, fill=PA)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod3")) +
  annotate("rect", xmax=15, xmin=6.608, ymin=-1, ymax=0.25, alpha=0.2,
           fill="black", color = "black", linetype = "dotted") +
  annotate("rect", xmax=-3.304, xmin=-15, ymin=-1, ymax=0.25, alpha=0.2,
           fill="black", color = "black", linetype = "dotted") +
  geom_vline(xintercept = outIQR2(lat.rs$asinhResilience), linetype = "longdash") +
  xlab("Recovery rate") +
  ylab("Density") + 
  coord_cartesian(ylim = c(0, 0.18), xlim = c(-10, 10)) +
  labs(fill = "Restoration\ntype") + 
  theme(legend.position = c(0.12,0.68), plot.margin = margins)
p3

#' # 3. Figure 1
#' 
#' Combining above plots into one figure and saving as .tiff (see R file for code)
#+ echo=FALSE
gp1a <- ggplotGrob(p1a)
gp1b <- ggplotGrob(p1b)
maxWidth = grid::unit.pmax(gp1a$widths[2:5], gp1b$widths[2:5])
gp1a$widths[2:5] <- as.list(maxWidth)
gp1b$widths[2:5] <- as.list(maxWidth)
gp2 <- ggplotGrob(p2)
gp3 <- ggplotGrob(p3)
maxWidth = grid::unit.pmax(gp2$widths[2:5], gp3$widths[2:5])
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

tiff(file = "Fig1.tif", width = 6, height = 10, units = "in", res = 600, 
     compression = "lzw")
grid.arrange(grid.arrange(grid.arrange(gp1a, gp1b, nrow = 2, heights = c(5,2)), 
                          gp2, gp3, nrow = 3, heights = c(12,7,7)))
grid.text(c("a.", "b.", "c."), x = 0.015, y = c(0.98, 0.53, 0.26), 
          vjust=1, hjust=0, gp=gpar(fontface=4,fontsize=13))
dev.off()

#' # 4. Statistical models
#' 
#' Results of all of the following models and sensitivity analyses are summarized 
#' in Table S1.
#' 
#' ## 4.1 Variance
#' 
#' JEA used used the rma.mv function (metafor) to analyze response metrics in a 
#' mixed-effects model framework (Viechtbauer 2010). rma.mv is designed to account 
#' for known, heteroscedastic sampling variances across studies and allows 
#' weighting by inverse variance, both of which are best-practices in meta-analysis 
#' (Gurevitch et al. 2018). However, JEA were unable to obtain the necessary variance 
#' data. They sought to compensate for this by setting variance to 1 for all data. 
#' Using a model designed for known, heteroscedastic variance by setting unknown 
#' variance to a uniform value violates a core assumption of the model. 
#' 
#' Given the absence of variance data, we considered it most appropriate to use 
#' generalized estimating equations (GEEs), which estimate population-averaged 
#' effects, do not require specification of observation-level variance, and produce 
#' standard errors robust to heteroscedasticity (Fieberg et al. 2009).
#' 
#' Another reasonable approach would be to use linear mixed-effects models (LMMs), 
#' which are broadly similar to the rma.mv model JEA used but entail estimation of 
#' unknown, observation-level variance rather than specification of known variance.
#' 
#' It is also possible that JEA's treatment of variance in rma.mv models, while 
#' inconsistent w/ model assumptions, does not qualitatively change interpretation.
#' 
#' Thus, to evaluate robustness of conclusions and their sensitivity to modeling 
#' decisions, we repeated all analyses using GEEs, LMMs, and rma.mv (w/ uniform 
#' variance following JEA).
#' 
#' - Fieberg J, RH Rieger, MC Zicus, JS Schildcrout. 2009 Regression modelling of correlated data in ecology: subject-specific and population averaged response patterns. J. Appl. Ecol. 46, 1018-1025.
#' - Gurevitch J, J Koricheva, S Nakagawa, G Stewart. 2018 Meta-analysis and the science of research synthesis. Nature 555, 175. (doi:10.1038/nature25753).
#' - Viechtbauer W. 2010 Conducting meta-analyses in R with the metafor package. 2010 36, 48. (doi:10.18637/jss.v036.i03).
#' 
#' ## 4.2 Random-error intercept
#' 
#' JEA included latitude as a random effect in addition to citation. But latitude 
#' was nearly completely confounded with citation. 
#' 
#' Number of distinct latitudes per citation
lat_cits <-
  latdata %>%
  group_by(Citation) %>%
  distinct(absLat) %>%
  summarise(n = n())
summary(lat_cits$n)

#' Proportion of citations w/ >1 latitude
lat_cits.sum <-
  lat_cits %>%
  filter(n > 1)
dim(lat_cits.sum)[1]/dim(lat_cits)[1]

#' Because of this overlap/redundancy between citation and latitude, we simplified 
#' our clusters/random-effects intercepts to include citation only. However, for 
#' comparability between JEA's rma.mv results and ours, we repeated all rma.mv 
#' models using both approaches (citation only and citation + latitude). 
#' 
#' ## 4.3 Outliers
#' 
#' While it is not uncommon to remove outliers in meta-analysis, it is critical to 
#' assess sensitivity of conclusions to outlier removal (Koricheva et al. 2013), 
#' which was not reported by JEA.
#' 
#' An outlier is defined as "an observation that lies outside the overall pattern 
#' of a distribution" (Moore and McCabe 2003). Examining JEA's outlier removal 
#' thresholds (as implemented and as intended), it does not appear that 
#' the discarded data were outside the overall distribution (figure 1c).
#'  
#' To assess sensitivity to outlier removal, we re-ran all models using complete 
#' data and using data with putative outliers removed following JEA's intended 
#' method. For comparison with JEA's original results, we additionally
#' implemented rma.mv models using the errant outlier removal thresholds (with 
#' latitude included as an additional random-error intercept). 
#' 
#' - Koricheva J, J Gurevitch, K Mengersen. 2013 Handbook of meta-analysis in ecology and evolution, Princeton University Press.
#' - Moore, D S, G P McCabe. 2003. Introduction to the Practice of Statistics. W.H. Freeman and Co. , New York.
#' 
#' ## 4.4 Generalized estimating equations
#' 
#' GEEs estimate average responses across populations rather than predicting effects 
#' at the observation (individual) level. Thus GEEs do not require knowledge of or 
#' estimation of observation-level variance. Standard errors from GEEs are robust 
#' to heteroscedasticity. 
#' 
#' Implemented using the "exchangeable" correlation structure---equivalent to 
#' assuming that observations from the same citation are equally correlated and 
#' that observations from different citations are independent (analogous to
#' including citation as random-error intercept in a mixed-effects model)
#' 
#' ### P-A only models (JEA "third model structure") 

#' Setting Passive as reference level in data so that all models show effect of 
#' Active restoration relative to Passive restoration 
lat.rs$PA <- factor(lat.rs$PA, levels(factor(lat.rs$PA))[c(2,1)])
lat.rs.or$PA <- factor(lat.rs.or$PA, levels(factor(lat.rs.or$PA))[c(2,1)])
lat.rs.or.manual$PA <- factor(lat.rs.or.manual$PA, levels(factor(lat.rs.or.manual$PA))[c(2,1)])

#' Complete data: PA highly significant
m.PA <- formula(asinhResilience ~ PA)
gee.PA <- geeglm(m.PA, data = lat.rs, id = Citation, corstr = "exchangeable")
(gee.PA.summary <- summary(gee.PA))
(gee.PA.anova <- anova(gee.PA))

#' Outlier removal: PA moderately significant
gee.PA.or <- geeglm(m.PA, data = lat.rs.or, id = Citation, corstr = "exchangeable")
(gee.PA.or.summary <- summary(gee.PA.or))
(gee.PA.or.anova <- anova(gee.PA.or))

#' ### PA + additional category models (JEA "first model structure") 
#' 
#' #### DisturbCat
#' 
#' Complete data: PA moderately significant, DisturbCat n.s.
m.PA.Dis <- formula(asinhResilience ~ DisturbCat + PA)
gee.PA.Dis <- geeglm(m.PA.Dis, data = lat.rs, id = Citation, corstr = "exchangeable")
(gee.PA.Dis.summary <- summary(gee.PA.Dis))
(gee.PA.Dis.anova <- anova(gee.PA.Dis))

#' Outlier removal: PA n.s., DisturbCat highly signficant
gee.PA.Dis.or <- geeglm(m.PA.Dis, data = lat.rs.or, id = Citation, corstr = "exchangeable")
(gee.PA.Dis.or.summary <- summary(gee.PA.Dis.or))
(gee.PA.Dis.or.anova <- anova(gee.PA.Dis.or))

#' #### HabitatCat
#' 
#' Complete data: PA and HabitatCat moderately significant
m.PA.Hab <- formula(asinhResilience ~ HabitatCat + PA)
gee.PA.Hab <- geeglm(m.PA.Hab, data = lat.rs, id = Citation, corstr = "exchangeable")
(gee.PA.Hab.summary <- summary(gee.PA.Hab))
(gee.PA.Hab.anova <- anova(gee.PA.Hab))

#' Outlier removal: PA n.s., HabitatCat highly signficant
gee.PA.Hab.or <- geeglm(m.PA.Hab, data = lat.rs.or, id = Citation, corstr = "exchangeable")
(gee.PA.Hab.or.summary <- summary(gee.PA.Hab.or))
(gee.PA.Hab.or.anova <- anova(gee.PA.Hab.or))

#' ## 4.5 Figure S1
#' 
#' ### Extracting estimates and CIs for figure
#' 
#' #### Third model structure, PA only (see R file for code)
#' All data
#+ echo=FALSE
newdat.PA <- expand.grid(PA = levels(lat.rs$PA), asinhResilience = 0)
mm.gee.PA = model.matrix(terms(gee.PA), newdat.PA)
newdat.PA$asinhResilience = mm.gee.PA %*% coef(gee.PA)
se.PA <- diag(mm.gee.PA %*% tcrossprod(gee.PA$geese$vbeta, mm.gee.PA))
newdat.PA <- data.frame(
  newdat.PA, Lower = newdat.PA$asinhResilience-2*sqrt(se.PA), 
  Upper = newdat.PA$asinhResilience+2*sqrt(se.PA),
  Factor = "Restoration type", Category = "Restoration type", Data = "All data")

#' Outliers removed
#+ echo=FALSE
newdat.PA.or <- expand.grid(PA = levels(lat.rs.or$PA), asinhResilience = 0)
mm.gee.PA.or = model.matrix(terms(gee.PA.or), newdat.PA)
newdat.PA.or$asinhResilience = mm.gee.PA.or %*% coef(gee.PA.or)
se.PA.or <- diag(mm.gee.PA.or %*% tcrossprod(gee.PA.or$geese$vbeta, mm.gee.PA.or))
newdat.PA.or <- data.frame(
  newdat.PA.or, Lower = newdat.PA.or$asinhResilience-2*sqrt(se.PA.or), 
  Upper = newdat.PA.or$asinhResilience+2*sqrt(se.PA.or),
  Factor = "Restoration type", Category = "Restoration type", Data = "Outlier removal")

#' #### First model structure, disturbance category (see R file for code)
#' All data 
#+ echo=FALSE
newdat.PA.Dis <- expand.grid(PA = levels(lat.rs$PA), DisturbCat = levels(lat.rs$DisturbCat), 
                             asinhResilience = 0)
mm.gee.PA.Dis = model.matrix(terms(gee.PA.Dis), newdat.PA.Dis)
newdat.PA.Dis$asinhResilience = mm.gee.PA.Dis %*% coef(gee.PA.Dis)
se.PA.Dis <- diag(mm.gee.PA.Dis %*% tcrossprod(gee.PA.Dis$geese$vbeta, mm.gee.PA.Dis))
newdat.PA.Dis <- data.frame(
  newdat.PA.Dis, Lower = newdat.PA.Dis$asinhResilience-2*sqrt(se.PA.Dis), 
  Upper = newdat.PA.Dis$asinhResilience+2*sqrt(se.PA.Dis), Factor = "Disturbance type", Data = "All data"
)
names(newdat.PA.Dis)[names(newdat.PA.Dis)=="DisturbCat"] <- "Category"

#' Outliers removed
#+ echo=FALSE
newdat.PA.Dis.or <- expand.grid(PA = levels(lat.rs.or$PA), DisturbCat = levels(lat.rs.or$DisturbCat), 
                                asinhResilience = 0)
mm.gee.PA.Dis.or = model.matrix(terms(gee.PA.Dis.or), newdat.PA.Dis.or)
newdat.PA.Dis.or$asinhResilience = mm.gee.PA.Dis.or %*% coef(gee.PA.Dis.or)
se.PA.Dis.or <- diag(mm.gee.PA.Dis.or %*% tcrossprod(gee.PA.Dis.or$geese$vbeta, mm.gee.PA.Dis.or))
newdat.PA.Dis.or <- data.frame(
  newdat.PA.Dis.or, Lower = newdat.PA.Dis.or$asinhResilience-2*sqrt(se.PA.Dis.or), 
  Upper = newdat.PA.Dis.or$asinhResilience+2*sqrt(se.PA.Dis.or), Factor = "Disturbance type", Data = "Outlier removal"
)
names(newdat.PA.Dis.or)[names(newdat.PA.Dis.or)=="DisturbCat"] <- "Category"

#' #### First model structure, habitat category (see R file for code)
#' All data 
#+ echo=FALSE
newdat.PA.Hab <- expand.grid(PA = levels(lat.rs$PA), HabitatCat = levels(lat.rs$HabitatCat), 
                             asinhResilience = 0)
mm.gee.PA.Hab = model.matrix(terms(gee.PA.Hab), newdat.PA.Hab)
newdat.PA.Hab$asinhResilience = mm.gee.PA.Hab %*% coef(gee.PA.Hab)
se.PA.Hab <- diag(mm.gee.PA.Hab %*% tcrossprod(gee.PA.Hab$geese$vbeta, mm.gee.PA.Hab))
newdat.PA.Hab <- data.frame(
  newdat.PA.Hab, Lower = newdat.PA.Hab$asinhResilience-2*sqrt(se.PA.Hab), 
  Upper = newdat.PA.Hab$asinhResilience+2*sqrt(se.PA.Hab), Factor = "Habitat type", Data = "All data"
)
names(newdat.PA.Hab)[names(newdat.PA.Hab)=="HabitatCat"] <- "Category"

#' Outliers removed 
#+ echo=FALSE
newdat.PA.Hab.or <- expand.grid(PA = levels(lat.rs.or$PA), HabitatCat = levels(lat.rs.or$HabitatCat), 
                                asinhResilience = 0)
mm.gee.PA.Hab.or = model.matrix(terms(gee.PA.Hab.or), newdat.PA.Hab.or)
newdat.PA.Hab.or$asinhResilience = mm.gee.PA.Hab.or %*% coef(gee.PA.Hab.or)
se.PA.Hab.or <- diag(mm.gee.PA.Hab.or %*% tcrossprod(gee.PA.Hab.or$geese$vbeta, mm.gee.PA.Hab.or))
newdat.PA.Hab.or <- data.frame(
  newdat.PA.Hab.or, Lower = newdat.PA.Hab.or$asinhResilience-2*sqrt(se.PA.Hab.or), 
  Upper = newdat.PA.Hab.or$asinhResilience+2*sqrt(se.PA.Hab.or), Factor = "Habitat type", Data = "Outlier removal"
)
names(newdat.PA.Hab.or)[names(newdat.PA.Hab.or)=="HabitatCat"] <- "Category"

#' Combining models into single table for plotting (see R file for code)
#+ echo=FALSE
gee.dat <- rbind(newdat.PA, newdat.PA.or, newdat.PA.Dis, newdat.PA.Dis.or, newdat.PA.Hab, newdat.PA.Hab.or)
gee.dat$asinhResilience <- as.vector(gee.dat$asinhResilience)

gee.ord <-
  gee.dat %>%
  filter(Data == "All data") %>%
  group_by(Factor, Category) %>%
  summarise(mean = mean(asinhResilience))  %>%
  arrange(mean)
gee.dat$Category <- factor(gee.dat$Category, as.character(gee.ord$Category))
gee.dat$PA <- factor(gee.dat$PA, levels(factor(gee.dat$PA))[c(2,1)])

#' ### Figure S1 (see R file for code)
#+ echo=FALSE
f5_labels <- data.frame(Data = rep(c("All data", "Outlier removal"), 2), 
                        label = c("'Active-Passive:'~italic('p')~'= 0.04'", 
                                  "'Active-Passive:'~italic('p')~'= 0.32'",
                                  "'Disturbance:'~italic('p')~'= 0.15'", 
                                  "'Disturbance:'~italic('p')~'= 0.002'"))
p5 <- ggplot(gee.dat[gee.dat$Factor=="Disturbance type",], 
             aes(x = Category, y = asinhResilience, color = PA)) +
  geom_point(aes(color = PA), size = 2, shape = 19, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), position = position_dodge(0.6), width = 0) +
  coord_flip() +
  scale_color_manual(values = c("dodgerblue4", "red4")) +
  xlab("Disturbance") + 
  ylab("Effect on recovery rate") +
  geom_text(y=2.5, x=c(2.3,1.7,2.3,1.7), aes(label=label), data=f5_labels, color = "black", hjust = 0, 
            parse = TRUE, size = 3) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  facet_grid(Data ~ Factor, scales = "free_y")

f6_labels <- data.frame(Data = rep(c("All data", "Outlier removal"), 2), 
                        label = c("'Active-Passive:'~italic('p')~'= 0.03'", 
                                  "'Active-Passive:'~italic('p')~'= 0.19'",
                                  "'Habitat:'~italic('p')~'= 0.02'", 
                                  "'Habitat:'~italic('p')~'= 0.0002'"))
p6 <- ggplot(gee.dat[gee.dat$Factor=="Habitat type",], aes(x = Category, y = asinhResilience, color = PA)) +
  geom_point(aes(color = PA), size = 2, shape = 19, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), position = position_dodge(0.6), width = 0) +
  coord_flip() +
  scale_color_manual(values = c("dodgerblue4", "red4")) +
  xlab("Habitat") + 
  ylab("Effect on recovery rate") +
  geom_text(y=3, x=c(2.3,1.7,2.3,1.7), aes(label=label), data=f6_labels, color = "black", hjust = 0, 
            parse = TRUE, size = 3) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  facet_grid(Data ~ Factor, scales = "free_y")

f7_labels <- data.frame(Data = c("All data", "Outlier removal"), 
                        label = c("italic('p')~'= 0.003'", "italic('p')~'= 0.04'"))
p7 <- ggplot(gee.dat[gee.dat$Factor=="Restoration type",], aes(x = Category, y = asinhResilience, color = PA)) +
  geom_point(aes(color = PA), size = 2, shape = 19, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), position = position_dodge(0.6), width = 0) +
  coord_flip() +
  scale_color_manual(values = c("dodgerblue4", "red4")) +
  ylab("Effect on recovery rate") +
  geom_text(y=1.05, x=1.35, aes(label=label), data=f7_labels, color = "black", hjust = 0, 
            parse = TRUE, size = 3) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.title.y=element_blank(), legend.position = "left") +
  labs(color="Rest.\ntype") +
  facet_wrap(Data ~ .)

gp5 <- ggplotGrob(p5)
gp6 <- ggplotGrob(p6)
gp7 <- ggplotGrob(p7)
maxWidth = grid::unit.pmax(gp5$widths[2:5])
gp5$widths[2:5] <- as.list(maxWidth)
gp6$widths[2:5] <- as.list(maxWidth)
#+ fig.height = 10
#+ echo=FALSE
grid.arrange(gp7, gp5, gp6, ncol=1, heights = c(1,3,3))

#+ echo=FALSE
tiff(file = "FigS1.tif", width = 6, height = 10, units = "in", res = 600, compression = "lzw")
grid.arrange(gp7, gp5, gp6, ncol=1, heights = c(1,3,3))
grid.text(c("a.","b.", "c."), x=0.125, y=c(0.99, 0.84, 0.42), vjust=1, hjust=0, gp=gpar(fontface=4))
dev.off()

#' ## 4.6 Linear mixed-effects models
#' 
#' LMMs estimate an additional variance parameter to account for variability at 
#' the observation level
#' 
#' ### PA-only models (JEA third model structure)
#' 
#' Complete data: PA highly significant
lmm.PA <- lmer(asinhResilience ~ PA + (1 | Citation), data = lat.rs)
(lmm.PA.summary <- summary(lmm.PA))
(lmm.PA.Anova <- Anova(lmm.PA))

#' Outlier removal: PA marginally significant
lmm.PA.or <- lmer(asinhResilience ~ PA + (1 | Citation), data = lat.rs.or)
(lmm.PA.or.summary <- summary(lmm.PA.or))
(lmm.PA.or.Anova <- Anova(lmm.PA.or))

#' ### PA + additional category models (JEA first model structure) 
#' 
#' #### DisturbCat
#'  
#' Complete data: PA moderately significant, DisturbCat n.s.
lmm.PA.Dis <- lmer(asinhResilience ~ PA + DisturbCat + (1 | Citation), data = lat.rs)
(lmm.PA.Dis.summary <- summary(lmm.PA.Dis))
(lmm.PA.Dis.Anova <- Anova(lmm.PA.Dis))

#' Outlier removal: PA n.s., DisturbCat moderately significant
lmm.PA.Dis.or <- lmer(asinhResilience ~ PA + DisturbCat + (1 | Citation), data = lat.rs.or)
lmm.PA.Dis.or
(lmm.PA.Dis.or.summary <- summary(lmm.PA.Dis.or))
(lmm.PA.Dis.or.Anova <- Anova(lmm.PA.Dis.or))

#' #### HabitatCat

#' Complete data: PA moderately significant, HabitatCat n.s.
lmm.PA.Hab <- lmer(asinhResilience ~ PA + HabitatCat + (1 | Citation), data = lat.rs)
(lmm.PA.Hab.summary <- summary(lmm.PA.Hab))
(lmm.PA.Hab.Anova <- Anova(lmm.PA.Hab))


#' Outlier removal: PA n.s., HabitatCat highly significant
lmm.PA.Hab.or <- lmer(asinhResilience ~ PA + HabitatCat + (1 | Citation), data = lat.rs.or)
(lmm.PA.Hab.or.summary <- summary(lmm.PA.Hab.or))
(lmm.PA.Hab.or.Anova <- Anova(lmm.PA.Hab.or))

#' ## 4.7 rma.mv models
#' 
#' Variance uniformly set to 1 following JEA. 
#' 
#' Using caching ('cache=TRUE') due to very long run times for these models.
#' 
#' ### Using simplified random-error intercept as above (Citation only)
#' 
#' ### PA-only models (JEA third model structure)
#' 
#' Complete data: PA highly significant
#+ cache=TRUE
Data = lat.rs
Data$vi = 1
rma.PA <- rma.mv(asinhResilience, vi, mods=~PA, random = list(~1|Citation), data = Data, method = "ML")
rma.PA

#' Outlier removal: PA moderately significant
#+ cache=TRUE
Data = lat.rs.or
Data$vi = 1
rma.PA.or <- rma.mv(asinhResilience, vi, mods=~PA, random = list(~1|Citation), data = Data, method = "ML")
rma.PA.or

#' ### PA + additional category models (JEA first model structure) 
#' 
#' #### DisturbCat
#' 
#' Complete data: PA highly significant, DisturbCat n.s.
#+ cache=TRUE
Data = lat.rs
Data$vi = 1
rma.PA.Dis <- rma.mv(asinhResilience, vi, mods=~DisturbCat + PA, random = list(~1|Citation), data = Data, method = "ML")
rma.PA.Dis
(rma.PA.Dis.anova <- anova(rma.PA, rma.PA.Dis)) #LRT for DisturbCat main effect

#' Outlier removal: PA moderately significant, DisturbCat highly significant
#+ cache=TRUE
Data = lat.rs.or
Data$vi = 1
rma.PA.Dis.or <- rma.mv(asinhResilience, vi, mods=~DisturbCat + PA, random = list(~1|Citation), data = Data, method = "ML")
rma.PA.Dis.or
(rma.PA.Dis.or.anova <- anova(rma.PA.or, rma.PA.Dis.or) )

#' #### HabitatCat
#'
#' Complete data: PA highly significant, HabitatCat marginally significant
#+ cache=TRUE
Data = lat.rs
Data$vi = 1
rma.PA.Hab <- rma.mv(asinhResilience, vi, mods=~HabitatCat + PA, random = list(~1|Citation), data = Data, method = "ML")
rma.PA.Hab
(rma.PA.Hab.anova <- anova(rma.PA, rma.PA.Hab) )

#' Outlier removal: PA moderately significant, DisturbCat highly significant
#+ cache=TRUE
Data = lat.rs.or
Data$vi = 1
rma.PA.Hab.or <- rma.mv(asinhResilience, vi, mods=~HabitatCat + PA, random = list(~1|Citation), data = Data, method = "ML")
rma.PA.Hab.or
(rma.PA.Hab.or.anova <- anova(rma.PA.or, rma.PA.Hab.or))

#' ### Using two-factor random-error intercept following JEA (Citation + Latitude)
#' 
#' ### PA-only models (JEA third model structure)
#' 
#' Complete data: PA highly significant
#+ cache=TRUE
Data = lat.rs
Data$vi = 1
rma.Lat.PA <- rma.mv(asinhResilience, vi, mods=~PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA

#' Outlier removal: PA moderately significant
#+ cache=TRUE
Data = lat.rs.or
Data$vi = 1
rma.Lat.PA.or <- rma.mv(asinhResilience, vi, mods=~PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.or

#' ### PA + additional category models (JEA first model structure) 
#' 
#' #### DisturbCat
#' 
#' Complete data: PA highly significant, DisturbCat n.s.
#+ cache=TRUE
Data = lat.rs
Data$vi = 1
rma.Lat.PA.Dis <- rma.mv(asinhResilience, vi, mods=~DisturbCat + PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.Dis
(rma.Lat.PA.Dis.anova <- anova(rma.Lat.PA, rma.Lat.PA.Dis)) #LRT for DisturbCat main effect

#' Outlier removal: PA moderately significant, DisturbCat highly significant
#+ cache=TRUE
Data = lat.rs.or
Data$vi = 1
rma.Lat.PA.Dis.or <- rma.mv(asinhResilience, vi, mods=~DisturbCat + PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.Dis.or
(rma.Lat.PA.Dis.or.anova <- anova(rma.Lat.PA.or, rma.Lat.PA.Dis.or) )

#' #### HabitatCat
#'
#' Complete data: PA highly significant, HabitatCat marginally significant
#+ cache=TRUE
Data = lat.rs
Data$vi = 1
rma.Lat.PA.Hab <- rma.mv(asinhResilience, vi, mods=~HabitatCat + PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.Hab
(rma.Lat.PA.Hab.anova <- anova(rma.Lat.PA, rma.Lat.PA.Hab) )

#' Outlier removal: PA moderately significant, DisturbCat highly significant
#+ cache=TRUE
Data = lat.rs.or
Data$vi = 1
rma.Lat.PA.Hab.or <- rma.mv(asinhResilience, vi, mods=~HabitatCat + PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.Hab.or
(rma.Lat.PA.Hab.or.anova <- anova(rma.Lat.PA.or, rma.Lat.PA.Hab.or))

#' ### Using errant outlier removal and two-factor random-error intercept following JEA
#' 
#' ### PA-only models (JEA third model structure)
#' 
#' Complete data: PA highly significant
#+ cache=TRUE

#' Errant outlier removal: 
#+ cache=TRUE
Data = lat.rs.or.manual
Data$vi = 1
rma.Lat.PA.or.err <- rma.mv(asinhResilience, vi, mods=~PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.or.err

#' ### PA + additional category models (JEA first model structure) 
#' 
#' #### DisturbCat
#' 
#' Errant outlier removal
#+ cache=TRUE
Data = lat.rs.or.manual
Data$vi = 1
rma.Lat.PA.Dis.or.err <- rma.mv(asinhResilience, vi, mods=~DisturbCat + PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.Dis.or.err
(rma.Lat.PA.Dis.or.err.anova <- anova(rma.Lat.PA.or.err, rma.Lat.PA.Dis.or.err) )

#' #### HabitatCat
#' 
#' Errant outlier removal
#+ cache=TRUE
Data = lat.rs.or.manual
Data$vi = 1
rma.Lat.PA.Hab.or.err <- rma.mv(asinhResilience, vi, mods=~HabitatCat + PA, random = list(~1|absLat, ~1|Citation), data = Data, method = "ML")
rma.Lat.PA.Hab.or.err
(rma.Lat.PA.Hab.or.err.anova <- anova(rma.Lat.PA.or.err, rma.Lat.PA.Hab.or.err))

#' ## 4.8 Table S1

#' Organizing GEEs (see R file for code)
#+ echo=FALSE
GEE.df <- rbind(data.frame(Model = "JEA_3rd", Data = "All", Method = "GEE", 
                           Variable = "ActiveRest", gee.PA.summary$coefficients[2,-3]),
                data.frame(Model = "JEA_3rd", Data = "IQR", Method = "GEE", 
                           Variable = "ActiveRest", gee.PA.or.summary$coefficients[2,-3]),
                data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "GEE", 
                           Variable = "DisturbCat", Estimate = NA, Std.err = NA,
                           "Pr(>|W|)" = gee.PA.Dis.anova$'P(>|Chi|)'[1]),
                data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "GEE", 
                           Variable = "ActiveRest", gee.PA.Dis.summary$coefficients
                           [row.names(gee.PA.Dis.summary$coefficients)=="PAActive",-3]),
                data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "GEE", 
                           Variable = "DisturbCat", Estimate = NA, Std.err = NA,
                           "Pr(>|W|)"  = gee.PA.Dis.or.anova$'P(>|Chi|)'[1]),
                data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "GEE", 
                           Variable = "ActiveRest", gee.PA.Dis.or.summary$coefficients
                           [row.names(gee.PA.Dis.or.summary$coefficients)=="PAActive",-3]),
                data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "GEE", 
                           Variable = "HabitatCat", Estimate = NA, Std.err = NA,
                           "Pr(>|W|)"  = gee.PA.Hab.anova$'P(>|Chi|)'[1]),
                data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "GEE", 
                           Variable = "ActiveRest", gee.PA.Hab.summary$coefficients
                           [row.names(gee.PA.Hab.summary$coefficients)=="PAActive",-3]),
                data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "GEE", 
                           Variable = "HabitatCat", Estimate = NA, Std.err = NA,
                           "Pr(>|W|)"  = gee.PA.Hab.or.anova$'P(>|Chi|)'[1]),
                data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "GEE", 
                           Variable = "ActiveRest", gee.PA.Hab.or.summary$coefficients
                           [row.names(gee.PA.Hab.or.summary$coefficients)=="PAActive",-3]))
colnames(GEE.df)[6:7] <- c("SE", "P")

#' Organizing LMMs (see R file for code)
#+ echo=FALSE
LMM.df <- rbind(data.frame(Model = "JEA_3rd", Data = "All", Method = "LMM", 
                           Variable = "ActiveRest",
                           Estimate = lmm.PA.summary$coefficients[2,'Estimate'],
                           SE = lmm.PA.summary$coefficients[2,'Std. Error'],
                           P = lmm.PA.Anova$'Pr(>Chisq)'),
                data.frame(Model = "JEA_3rd", Data = "IQR", Method = "LMM", 
                           Variable = "ActiveRest",
                           Estimate = lmm.PA.or.summary$coefficients[2,'Estimate'],
                           SE = lmm.PA.or.summary$coefficients[2,'Std. Error'],
                           P = lmm.PA.or.Anova$'Pr(>Chisq)'),
                data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "LMM", 
                           Variable = "DisturbCat", Estimate = NA, SE = NA,
                           P = lmm.PA.Dis.Anova$'Pr(>Chisq)'[2]),
                data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "LMM", 
                           Variable = "ActiveRest", Estimate = lmm.PA.Dis.summary$coefficients
                           [row.names(lmm.PA.Dis.summary$coefficients)=="PAActive",'Estimate'],
                           SE = lmm.PA.Dis.summary$coefficients
                           [row.names(lmm.PA.Dis.summary$coefficients)=="PAActive",'Std. Error'],
                           P = lmm.PA.Dis.Anova$'Pr(>Chisq)'[1]),
                data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "LMM", 
                           Variable = "DisturbCat", Estimate = NA, SE = NA, 
                           P = lmm.PA.Dis.or.Anova$'Pr(>Chisq)'[2]),
                data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "LMM", 
                           Variable = "ActiveRest", Estimate = lmm.PA.Dis.or.summary$coefficients
                           [row.names(lmm.PA.Dis.or.summary$coefficients)=="PAActive",'Estimate'],
                           SE = lmm.PA.Dis.or.summary$coefficients
                           [row.names(lmm.PA.Dis.or.summary$coefficients)=="PAActive",
                             'Std. Error'],
                           P = lmm.PA.Dis.or.Anova$'Pr(>Chisq)'[1]),
                data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "LMM", 
                           Variable = "HabitatCat", Estimate = NA, SE = NA, 
                           P = lmm.PA.Hab.Anova$'Pr(>Chisq)'[2]),
                data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "LMM", 
                           Variable = "ActiveRest", Estimate = lmm.PA.Hab.summary$coefficients
                           [row.names(lmm.PA.Hab.summary$coefficients)=="PAActive",'Estimate'],
                           SE = lmm.PA.Hab.summary$coefficients
                           [row.names(lmm.PA.Hab.summary$coefficients)=="PAActive",'Std. Error'],
                           P = lmm.PA.Hab.Anova$'Pr(>Chisq)'[1]),
                data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "LMM", 
                           Variable = "HabitatCat", Estimate = NA, SE = NA,
                           P = lmm.PA.Hab.or.Anova$'Pr(>Chisq)'[2]),
                data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "LMM", 
                           Variable = "ActiveRest", Estimate = lmm.PA.Hab.or.summary$coefficients
                           [row.names(lmm.PA.Hab.or.summary$coefficients)=="PAActive",'Estimate'],
                           SE = lmm.PA.Hab.or.summary$coefficients
                           [row.names(lmm.PA.Hab.or.summary$coefficients)=="PAActive",
                             'Std. Error'],
                           P = lmm.PA.Hab.or.Anova$'Pr(>Chisq)'[1]))

#' Organizing rma.mv (see R file for code)
#+ echo=FALSE
RMA.df <- rbind(data.frame(Model = "JEA_3rd", Data = "All", Method = "RMA", 
                           Variable = "ActiveRest", data.frame(rma.PA[c(2,3,5)])[-1,]),
                data.frame(Model = "JEA_3rd", Data = "IQR", Method = "RMA", 
                           Variable = "ActiveRest", data.frame(rma.PA.or[c(2,3,5)])[-1,]),
                data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "RMA", 
                           Variable = "ActiveRest", rma.PA.Dis[c(2,3,5)])["PAActive",],
                data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "RMA", 
                           Variable = "DisturbCat", beta = NA, se = NA, 
                           pval = rma.PA.Dis.anova$pval),
                data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "RMA", 
                           Variable = "ActiveRest", rma.PA.Dis.or[c(2,3,5)])["PAActive",],
                data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "RMA", 
                           Variable = "DisturbCat", beta = NA, se = NA, 
                           pval = rma.PA.Dis.or.anova$pval),
                data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "RMA", 
                           Variable = "ActiveRest", rma.PA.Hab[c(2,3,5)])["PAActive",],
                data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "RMA", 
                           Variable = "HabitatCat", beta = NA, se = NA, 
                           pval = rma.PA.Hab.anova$pval),
                data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "RMA", 
                           Variable = "ActiveRest", rma.PA.Hab.or[c(2,3,5)])["PAActive",],
                data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "RMA", 
                           Variable = "HabitatCat", beta = NA, se = NA, 
                           pval = rma.PA.Hab.or.anova$pval))
colnames(RMA.df)[5:7] <- c("Estimate", "SE", "P")

RMA.Lat.df <- rbind(data.frame(Model = "JEA_3rd", Data = "All", Method = "RMA.Lat", 
                               Variable = "ActiveRest", data.frame(rma.Lat.PA[c(2,3,5)])[-1,]),
                    data.frame(Model = "JEA_3rd", Data = "IQR", Method = "RMA.Lat", 
                               Variable = "ActiveRest", data.frame(rma.Lat.PA.or[c(2,3,5)])[-1,]),
                    data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "RMA.Lat", 
                               Variable = "ActiveRest", rma.Lat.PA.Dis[c(2,3,5)])["PAActive",],
                    data.frame(Model = "JEA_1st_Dis", Data = "All", Method = "RMA.Lat", 
                               Variable = "DisturbCat", beta = NA, se = NA, 
                               pval = rma.Lat.PA.Dis.anova$pval),
                    data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "RMA.Lat", 
                               Variable = "ActiveRest", rma.Lat.PA.Dis.or[c(2,3,5)])["PAActive",],
                    data.frame(Model = "JEA_1st_Dis", Data = "IQR", Method = "RMA.Lat", 
                               Variable = "DisturbCat", beta = NA, se = NA, 
                               pval = rma.Lat.PA.Dis.or.anova$pval),
                    data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "RMA.Lat", 
                               Variable = "ActiveRest", rma.Lat.PA.Hab[c(2,3,5)])["PAActive",],
                    data.frame(Model = "JEA_1st_Hab", Data = "All", Method = "RMA.Lat", 
                               Variable = "HabitatCat", beta = NA, se = NA, 
                               pval = rma.Lat.PA.Hab.anova$pval),
                    data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "RMA.Lat", 
                               Variable = "ActiveRest", rma.Lat.PA.Hab.or[c(2,3,5)])["PAActive",],
                    data.frame(Model = "JEA_1st_Hab", Data = "IQR", Method = "RMA.Lat", 
                               Variable = "HabitatCat", beta = NA, se = NA, 
                               pval = rma.Lat.PA.Hab.or.anova$pval))
colnames(RMA.Lat.df)[5:7] <- c("Estimate", "SE", "P")

RMA.Lat.err.df <- rbind(data.frame(Model = "JEA_3rd", Data = "Errant.OR", Method = "RMA.Lat", 
                                   Variable = "ActiveRest", 
                                   data.frame(rma.Lat.PA.or.err[c(2,3,5)])[-1,]),
                        data.frame(Model = "JEA_1st_Dis", 
                                   Data = "Errant.OR", Method = "RMA.Lat", 
                                   Variable = "ActiveRest", rma.Lat.PA.Dis.or.err[c(2,3,5)])["PAActive",],
                        data.frame(Model = "JEA_1st_Dis", 
                                   Data = "Errant.OR", Method = "RMA.Lat", 
                                   Variable = "DisturbCat", beta = NA, se = NA, 
                                   pval = rma.Lat.PA.Dis.or.err.anova$pval),
                        data.frame(Model = "JEA_1st_Hab", 
                                   Data = "Errant.OR", Method = "RMA.Lat", 
                                   Variable = "ActiveRest", rma.Lat.PA.Hab.or.err[c(2,3,5)])["PAActive",],
                        data.frame(Model = "JEA_1st_Hab", 
                                   Data = "Errant.OR", Method = "RMA.Lat", 
                                   Variable = "HabitatCat", beta = NA, se = NA, 
                                   pval = rma.Lat.PA.Hab.or.err.anova$pval))
colnames(RMA.Lat.err.df)[5:7] <- c("Estimate", "SE", "P")

#' Combining all models
Models.df <- rbind(GEE.df, LMM.df, RMA.df, RMA.Lat.df, RMA.Lat.err.df)
Models.df <- arrange(Models.df, Model, Data, Method, Variable)
Models.df$Significance <- 0
Models.df$Significance[Models.df$P > 0.1] <- ""
Models.df$Significance[Models.df$P <= 0.1] <- "."
Models.df$Significance[Models.df$P <= 0.05] <- "*"
Models.df$Significance[Models.df$P <= 0.01] <- "**"
Models.df$Significance[Models.df$P <= 0.001] <- "***"
Models.df

write.csv(Models.df, file = "Models.df.csv", row.names = FALSE)

#' ### R Information (version and packages)
devtools::session_info()
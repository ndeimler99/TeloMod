#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

AXIS_TITLE_SIZE <- 11
AXIS_TEXT_SIZE <- 9
# 
# args <- c("/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/TELOMERE_MODS/131.5mC.summary.txt",
#           "/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/131_TARPON/sample/sample.telo_stats.txt",
#           "/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/CLUSTERING/1000.2000.cluster_assignment.txt")

args <- c("/media/data_01/ndeimler/HUMAN/TELOMOD/131_HEK/work/93/12611b4e9134180c8a5cd32e25df0d/telomeric_mod_summary.txt",
          "/media/data_01/ndeimler/HUMAN/TELOMOD/131_HEK/work/93/12611b4e9134180c8a5cd32e25df0d/sample.telo_stats.txt", FALSE)
#           "/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/CLUSTERING/1000.2000.cluster_assignment.txt")

args <- commandArgs(trailing=TRUE)
flag <- tolower(args[3]) == "true"
##### Plot Clustering Results #####

df <- read.table(args[2], header=TRUE)
df <- df %>% group_by(Cluster) %>% summarize(count=n())

plt <- ggplot(data=df) +
  geom_boxplot(mapping=aes(x=1, y=count/sum(count) * 100)) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("Telomeric Reads per Cluster [%]") +
  geom_hline(mapping=aes(yintercept=1.08), color="red", linetype="dashed")
plt
ggsave("cluster_size_distribution.pdf", plot=plt, width=12, height=10, units="cm")

#####  Telomere Modification Stats in Cluster Specific Manner #####
telo_stats <- read.table(args[2], header=TRUE)
telo_modifications <- read.table(args[1], header=TRUE)
telo_modifications <- telo_modifications %>% dplyr::select("read_id", "read_type", "proportion")
telo_modifications <- telo_modifications %>% pivot_wider(names_from = read_type, values_from = proportion)

telo_stats <- left_join(telo_stats, telo_modifications, by="read_id")

telo_stats$Cluster <- factor(telo_stats$Cluster)
telo_stats <- telo_stats[!is.na(telo_stats$Cluster),]

order_df <- telo_stats %>% group_by(Cluster) %>% 
  summarise(median_len=median(vrr_telo_length), .groups="drop") %>%
  arrange(median_len) %>%
  mutate(cluster=factor(Cluster, levels=Cluster))
telo_stats$Cluster <- factor(telo_stats$Cluster, levels=levels(order_df$Cluster))


plt <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=Cluster, y=subtelo), linewidth=0.3, outlier.size=1) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Cluster") + ylab("Subtelomere Modifications [%]")

ggsave("cluster_subtelo_modifications.pdf", plot=plt, width=25, height=10, units="cm")

if (!flag){
  plt <- ggplot(data=telo_stats) +
    geom_boxplot(mapping=aes(x=Cluster, y=subtelo, fill=strand), linewidth=0.3, outlier.size=1) +
    theme_minimal() +
    theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
          axis.text=element_text(size=AXIS_TEXT_SIZE),
          axis.text.x=element_text(angle=45, hjust=1)) +
    xlab("Cluster") + ylab("Subtelomere Modifications [%]") +
    scale_fill_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand")

  ggsave("cluster_by_strand_subtelo_modifications.pdf", plot=plt, width=25, height=10, units="cm")
}

plt <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=Cluster, y=telo), linewidth=0.3, outlier.size=1) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Cluster") + ylab("Subtelomere Modifications [%]")

ggsave("cluster_telo_modifications.pdf", plot=plt, width=25, height=10, units="cm")

if (!flag){
  plt <- ggplot(data=telo_stats) +
    geom_boxplot(mapping=aes(x=Cluster, y=telo, fill=strand), linewidth=0.3, outlier.size=1) +
    theme_minimal() +
    theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
          axis.text=element_text(size=AXIS_TEXT_SIZE),
          axis.text.x=element_text(angle=45, hjust=1)) +
    xlab("Cluster") + ylab("Subtelomere Modifications [%]") +
    scale_fill_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand")

  ggsave("cluster_by_strand_telo_modifications.pdf", plot=plt, width=25, height=10, units="cm")
}


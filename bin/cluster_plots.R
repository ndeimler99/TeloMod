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


args <- commandArgs(trailing=TRUE)

##### Plot Clustering Results #####

df <- read.table(args[3], header=TRUE)
df <- df %>% group_by(cluster) %>% summarize(count=n())

plt <- ggplot(data=df) +
  geom_boxplot(mapping=aes(x=1, y=count/sum(count) * 100)) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("Telomeric Reads per Cluster [%]") +
  geom_hline(mapping=aes(yintercept=1.08), color="red", linetype="dashed")

ggsave("cluster_size_distribution.pdf", plot=plt, width=6, height=4, units="cm")

#####  Telomere Modification Stats in Cluster Specific Manner #####
cluster_assignment <- read.table(args[3], header=TRUE)
telo_stats <- read.table(args[2], header=TRUE)
telo_modifications <- read.table(args[1], header=TRUE)
telo_modifications <- telo_modifications %>% dplyr::select("read_id", "read_type", "proportion")
telo_modifications <- telo_modifications %>% pivot_wider(names_from = read_type, values_from = proportion)

telo_stats <- left_join(telo_stats, cluster_assignment, by="read_id")
telo_stats <- left_join(telo_stats, telo_modifications, by="read_id")

telo_stats$cluster <- factor(telo_stats$cluster)
telo_stats <- telo_stats[!is.na(telo_stats$cluster),]

order_df <- telo_stats %>% group_by(cluster) %>% 
  summarise(median_len=median(vrr_telo_length), .groups="drop") %>%
  arrange(median_len) %>%
  mutate(cluster=factor(cluster, levels=cluster))
telo_stats$cluster <- factor(telo_stats$cluster, levels=levels(order_df$cluster))


plt <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=cluster, y=subtelo)) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Cluster") + ylab("Subtelomere Modifications [%]")
plt
ggsave("cluster_subtelo_modifications.pdf", plot=plt, width=12, height=5, units="cm")

plt <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=cluster, y=subtelo, fill=strand)) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Cluster") + ylab("Subtelomere Modifications [%]") +
  scale_fill_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand")

ggsave("cluster_by_strand_subtelo_modifications.pdf", plot=plt, width=12, height=5, units="cm")
  

plt <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=cluster, y=telo)) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Cluster") + ylab("Subtelomere Modifications [%]")

ggsave("cluster_telo_modifications.pdf", plot=plt, width=12, height=5, units="cm")

plt <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=cluster, y=telo, fill=strand)) +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Cluster") + ylab("Subtelomere Modifications [%]") +
  scale_fill_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand")

ggsave("cluster_by_strand_telo_modifications.pdf", plot=plt, width=12, height=5, units="cm")



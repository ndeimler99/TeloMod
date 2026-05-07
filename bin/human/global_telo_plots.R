#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

AXIS_TITLE_SIZE <- 11
AXIS_TEXT_SIZE <- 9

# args <- c("/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/131_HEK_5mC/genomic_summary.txt",
#           "/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/131_HEK_5mC/telomeric_mod_summary.txt",
#           "/media/data_01/ndeimler/TELOMERE_MODIFICATIONS/131_HEK/131_TARPON/sample/sample.telo_stats.txt")

args <- commandArgs(trailing=TRUE)

flag <- tolower(args[3]) == "true"

##### Load Relevant Data ####
telomeric_df <- read.table(args[1], header=TRUE)

telomeric_df <- telomeric_df %>% dplyr::select("read_id", "strand", "read_type", "proportion")
colnames(telomeric_df) <- c("read_id", "strand", "read_type", "mod_percentage")
df <- telomeric_df
df <- df[df$mod_percentage <=100,]
df$read_type <- factor(df$read_type, levels=c("full_read", "subtelo", "telo"))

##### modification comparison between Genomic and Telomeric Reads ####

plt <- ggplot(data=df[df$read_type %in% c("full_read"),]) +
  geom_boxplot(mapping=aes(x=read_type, y=mod_percentage), outliers=FALSE) +
  theme_minimal() +
  xlab("Read Type") + ylab("Nucleotides Modified [%]") +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE)) +
  scale_x_discrete(breaks=c("full_read"), labels=c("Telomeric"))

ggsave("telomeric.pdf", plot=plt, width=12, height=10, units="cm")

#####  modification comparison between Genomic, Telomeric Reads split by reasoning ##### 
plt <- ggplot(data=df[df$read_type != "full_read",]) +
  geom_boxplot(mapping=aes(x=read_type, y=mod_percentage), outliers=FALSE) +
  theme_minimal() +
  xlab("Read Type") + ylab("Nucleotides Modified [%]") +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE)) +
  scale_x_discrete(breaks=c("subtelo", "telo"), labels=c("Subtelo", "Telomere"))

ggsave("telomeric_segregated_by_type.pdf", plot=plt, width=12, height=10, units="cm")

if (!flag){
  #####  modification comparison between Genomic, Telomeric Reads split by reasoning and segregated by strand
  plt <- ggplot(data=df[df$read_type != "full_read",]) +
    geom_boxplot(mapping=aes(x=read_type, y=mod_percentage, fill=strand), outliers=FALSE) +
    theme_minimal() +
    xlab("Read Type") + ylab("Nucleotides Modified [%]") +
    theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
          axis.text=element_text(size=AXIS_TEXT_SIZE)) +
    scale_x_discrete(breaks=c("subtelo", "telo"), labels=c("Subtelo", "Telomere")) +
    scale_fill_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"))

  ggsave("telomeric_segregated_by_type_and_strand.pdf", plot=plt, width=12, height=10, units="cm")
}

#####  modification percentage by read length ##### 
telo_stats <- read.table(args[2], header=TRUE)
telomeric_df <- telomeric_df %>% pivot_wider(names_from = read_type, values_from = mod_percentage)
telo_stats <- left_join(telo_stats, telomeric_df, by="read_id")


plt <- ggplot(data=telo_stats) +
  geom_point(mapping=aes(x=telo, y=read_len, color=strand.x), alpha=0.5) +
  scale_color_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand") +
  xlab("Telomeric Modifications [%]") + ylab("Read Length") +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE))

ggsave("telo_mods_vs_read_length.pdf", plot=plt, width=12, height=10, units="cm")

plt <- ggplot(data=telo_stats) +
  geom_point(mapping=aes(x=subtelo, y=read_len, color=strand.x), alpha=0.5) +
  scale_color_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand") +
  xlab("Sub-Telomeric Modifications [%]") + ylab("Read Length") +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE))

ggsave("subtelo_mods_vs_read_length.pdf", plot=plt, width=12, height=10, units="cm")

#####  modification percentage by subtelomere length ##### 

plt <- ggplot(data=telo_stats) +
  geom_point(mapping=aes(x=telo, y=read_len-vrr_telo_length, color=strand.x), alpha=0.5) +
  scale_color_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand") +
  xlab("Telomeric Modifications [%]") + ylab("Subtelomere Length") +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE))
ggsave("telo_mods_vs_subtelo_length.pdf", plot=plt, width=12, height=10, units="cm")

plt <- ggplot(data=telo_stats) +
  geom_point(mapping=aes(x=subtelo, y=read_len-vrr_telo_length, color=strand.x), alpha=0.5) +
  scale_color_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand") +
  xlab("Sub-Telomeric Modifications [%]") + ylab("Subtelomere Length") +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE))
ggsave("subtelo_mods_vs_subtelo_length.pdf", plot=plt, width=12, height=10, units="cm")

#####  modification percentage by telomere length ##### 

plt <- ggplot(data=telo_stats) +
  geom_point(mapping=aes(x=telo, y=vrr_telo_length, color=strand.x), alpha=0.5) +
  scale_color_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand") +
  xlab("Telomeric Modifications [%]") + ylab("Telomere Length") +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE))
ggsave("telo_mods_vs_telo_length.pdf", plot=plt, width=12, height=10, units="cm")

plt <- ggplot(data=telo_stats) +
  geom_point(mapping=aes(x=subtelo, y=vrr_telo_length, color=strand.x), alpha=0.5) +
  scale_color_manual(breaks=c("C", "G"), values=c("#D81B60", "#1E88E5"), name="Strand") +
  xlab("Sub-Telomeric Modifications [%]") + ylab("Telomere Length") +
  theme_minimal() +
  theme(axis.title=element_text(size=AXIS_TITLE_SIZE),
        axis.text=element_text(size=AXIS_TEXT_SIZE))
ggsave("subtelo_mods_vs_telo_length.pdf", plot=plt, width=12, height=10, units="cm")



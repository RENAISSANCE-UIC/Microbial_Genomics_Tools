
# Load libraries
library(ggplot2)
library(dplyr)

# Read coverage data
setwd("/home/william-ackerman/Desktop/cnMOPS/read_depth")

list.files(pattern = "_depth.txt")

cov_all <- bind_rows(
  read.table("A1-D10_S11_depth.txt") %>% mutate(replicate = "A1-D10"),
  read.table("A2-D10_S12_depth.txt") %>% mutate(replicate = "A2-D10"),
  read.table("A3-D10_S7_depth.txt") %>% mutate(replicate = "A3-D10"),
  read.table("A4-D10_S15_depth.txt") %>% mutate(replicate = "A4-D10"),
  read.table("A5-D10_S13_depth.txt") %>% mutate(replicate = "A5-D10")
)
colnames(cov_all) <- c("chrom", "pos", "depth", "replicate")

# Normalize each replicate by its own genome-wide mean
cov_all <- cov_all %>%
  group_by(replicate) %>%
  mutate(norm_depth = depth / mean(depth))

# Plot
# ggplot(cov_all, aes(x = pos, y = norm_depth, color = replicate)) +
#   geom_line() +
#   geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray") +
#   geom_hline(yintercept = 1.8, linetype = "dashed", color = "red") +
#   labs(title = "Coverage comparison across replicates",
#        x = "Position", y = "Normalized coverage") +
#   theme_minimal()


# Faceted plot
ggplot(cov_all, aes(x = pos, y = norm_depth)) +
  geom_line(color = "steelblue4", alpha = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray43") +
  geom_hline(yintercept = 1.8, linetype = "dashed", color = "firebrick4") +
  facet_wrap(~ replicate, ncol = 2) +
  labs(title = "Normalized Coverage Across Region by Replicate",
       x = "Position", y = "Normalized Coverage") +
  ylim(c(0,5)) +
  DOSE::theme_dose() 


# Untreated (no selective pressure)
cov_all2 <- bind_rows(
  read.table("U1-D10_S6_depth.txt") %>% mutate(replicate = "U1-D10"),
  read.table("U2-D10_S7_depth.txt") %>% mutate(replicate = "U2-D10"),
  read.table("U3-D10_S8_depth.txt") %>% mutate(replicate = "U3-D10"),
  read.table("U4-D10_S9_depth.txt") %>% mutate(replicate = "U4-D10"),
  read.table("U5-D10_S10_depth.txt") %>% mutate(replicate = "U5-D10")
)
colnames(cov_all2) <- c("chrom", "pos", "depth", "replicate")

# Normalize each replicate by its own genome-wide mean
cov_all2 <- cov_all2 %>%
  group_by(replicate) %>%
  mutate(norm_depth = depth / mean(depth))

# Faceted plot
ggplot(cov_all2, aes(x = pos, y = norm_depth)) +
  geom_line(color = "steelblue4", alpha = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray43") +
  geom_hline(yintercept = 1.8, linetype = "dashed", color = "firebrick4") +
  facet_wrap(~ replicate, ncol = 2) +
  labs(title = "Normalized Coverage Across Region by Replicate",
       x = "Position", y = "Normalized Coverage") +
  ylim(c(0,5)) +
  DOSE::theme_dose() 


# Cephalexin treated 
cov_all3 <- bind_rows(
  read.table("C1-D10_S1_depth.txt") %>% mutate(replicate = "C1-D10"),
  read.table("C2-D10_S2_depth.txt") %>% mutate(replicate = "C2-D10"),
  read.table("C3-D10_S3_depth.txt") %>% mutate(replicate = "C3-D10"),
  read.table("C4-D10_S4_depth.txt") %>% mutate(replicate = "C4-D10"),
  read.table("C5-D10_S5_depth.txt") %>% mutate(replicate = "C5-D10")
)
colnames(cov_all3) <- c("chrom", "pos", "depth", "replicate")

# Normalize each replicate by its own genome-wide mean
cov_all3 <- cov_all3 %>%
  group_by(replicate) %>%
  mutate(norm_depth = depth / mean(depth))

# Faceted plot
ggplot(cov_all3, aes(x = pos, y = norm_depth)) +
  geom_line(color = "steelblue4", alpha = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray43") +
  geom_hline(yintercept = 1.8, linetype = "dashed", color = "firebrick4") +
  facet_wrap(~ replicate, ncol = 2) +
  labs(title = "Normalized Coverage Across Region by Replicate",
       x = "Position", y = "Normalized Coverage") +
  ylim(c(0,5)) +
  DOSE::theme_dose() 

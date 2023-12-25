#Libraries
library(ggplot2)
library(dplyr)
library(rstatix)


#Figure 1####
#Panel A

# Perform pairwise t-tests
pairwise_results <- fig1a_long %>%
  pairwise_t_test(Intensity ~ condition, p.adjust.method = "bonferroni", ref.group = "Un.")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig1A_stats.csv')

#Panel B
pairwise_results <- fig1b_long %>%
  pairwise_t_test(density ~ condition, p.adjust.method = "bonferroni", ref.group = "Un.")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig1B_stats.csv')


#Panel C
# Perform pairwise t-tests within each level of Condition
pairwise_results <- fig1c_long %>%
  group_by(Condition) %>%
  pairwise_t_test(Intensity ~ Treatment, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig1C_stats.csv')

#Panel D
pairwise_results <- fig1d_long %>%
  group_by(Condition) %>%
  pairwise_t_test(Intensity ~ Treatment, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig1D_stats.csv')

#Figure 2####
#Panel A
pairwise_results <- fig2a_long %>%
  pairwise_t_test(density ~ condition, p.adjust.method = "bonferroni", ref.group = "HaCaTWT")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig2A_stats.csv')

#Panel B
pairwise_results <- fig2b_long %>%
  pairwise_t_test(density ~ condition, p.adjust.method = "bonferroni", ref.group = "HaCaTWT")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig2B_stats.csv')

#Panel C
pairwise_results <- fig2c_long %>%
  group_by(Condition) %>% 
  pairwise_t_test(Intensity ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig2C_stats.csv')

#Panel D
pairwise_results <- fig2d_long %>%
  group_by(Condition) %>% 
  pairwise_t_test(Intensity ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig2D_stats.csv')


#Panel E
pairwise_results <- fig2e_long %>%
  group_by(Condition) %>% 
  pairwise_t_test(Intensity ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig2E_stats.csv')
#Panel F
pairwise_results <- fig2f_long %>%
  group_by(Condition) %>% 
  pairwise_t_test(Intensity ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig2F_stats.csv')
#Figure 3####
#Panel A
pairwise_results <- fig3a_long %>%
  group_by(Condition) %>% 
  pairwise_t_test(Intensity ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig3A_stats.csv')

#Panel B
pairwise_results <- fig3b_long %>%
  group_by(Condition) %>% 
  pairwise_t_test(Intensity ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig3B_stats.csv')
#Figure 4####
#Panel B
pairwise_results <- fig4b_long %>%
  group_by(condition) %>% 
  pairwise_t_test(count ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig4B_stats.csv')
#Panel C
pairwise_results <- fig4c_long %>%
  pairwise_t_test(area ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig4C_stats.csv')


#Panel D
pairwise_results <- fig4d_long %>%
  pairwise_t_test(area ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig4D_stats.csv')

#Panel E
pairwise_results <- fig4e_long %>%
  pairwise_t_test(area ~ cell_type, p.adjust.method = "bonferroni")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig4E_stats.csv')

#Figure 5####
#Panel A
pairwise_results <- fig5a_long %>%
  group_by(Time) %>% 
  pairwise_t_test(titre ~ cell_type, p.adjust.method = "bonferroni", ref.group = "HaCaTWT")

# View the results
print(pairwise_results)

write.csv(pairwise_results, 'Fig5A_stats.csv')


#Panel C
results <- pairwise.t.test(data_uniqueID_CP_metrics$TrackObjects_IntegratedDistance, 
                           data_uniqueID_CP_metrics$Metadata_Condition, 
                           p.adjust.method = "bonferroni")
print(results)

write.csv(pairwise_results, 'Fig5C_stats.csv')

#Panel D
results <- pairwise.t.test(directionality_df$directionality, 
                           directionality_df$Metadata_Condition.y, 
                           p.adjust.method = "bonferroni")
print(results)

write.csv(pairwise_results, 'Fig5D_stats.csv')

#Panel E
results <- pairwise.t.test(data_uniqueID_per_cell$radial_velocity, 
                           data_uniqueID_per_cell$Metadata_Condition, 
                           p.adjust.method = "bonferroni")
print(results)

write.csv(pairwise_results, 'Fig5E_stats.csv')

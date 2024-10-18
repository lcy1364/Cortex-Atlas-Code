library(tidyverse)
set("~/cortex/figS1-6/")
wholeMeta <- read.csv("wholeMeta.csv")

# IT
IT <- wholeMeta %>% filter(str_detect(cross_area_subclass, "IT")) %>%  group_by(donor, batch, cross_area_subclass) %>% summarise(count = n()) %>% group_by(donor, batch) %>% mutate(per = count /
                                                                                                                                                                                      sum(count)) %>% group_by(batch, cross_area_subclass) %>% summarise(aver = median(per)) %>% pivot_wider(names_from = batch, values_from = aver)

md <- lm(edlein ~ us, IT)

# > summary(md)
#
# Call:
#   lm(formula = edlein ~ us, data = IT)
#
# Residuals:
#   1         2         3         4         5
# -0.006503  0.020444  0.012081  0.011364 -0.037387
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.05435    0.01547   3.512   0.0391 *
#   us           0.71034    0.04873  14.577   0.0007 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.02667 on 3 degrees of freedom
# Multiple R-squared:  0.9861,	Adjusted R-squared:  0.9814
# F-statistic: 212.5 on 1 and 3 DF,  p-value: 0.0007001


p1 <- ggplot(IT, mapping = aes(x = edlein, y = us)) + geom_point(color = "grey30") + geom_smooth(method = "lm", se = T) + coord_fixed(xlim = c(0, 0.7),
                                                                                                                                      ylim = c(0, 0.7),
                                                                                                                                      expand = F) + cowplot::theme_cowplot() + labs(title = "subclass",subtitle = "Adjusted R-squared:  0.9814\n P = 0.0007")


# GABA
GABA <- wholeMeta %>% filter(str_detect(class, "inhibitory")) %>%  group_by(donor, batch, cross_area_subclass) %>% summarise(count = n()) %>% group_by(donor, batch) %>% mutate(per = count /
                                                                                                                                                                                  sum(count)) %>% group_by(batch, cross_area_subclass) %>% summarise(aver = median(per)) %>% pivot_wider(names_from = batch, values_from = aver)

md <- lm(edlein ~ us, GABA)

# summary(md)
# Call:
#   lm(formula = edlein ~ us, data = GABA)
#
# Residuals:
#   Min        1Q    Median        3Q       Max
# -0.043235 -0.014162 -0.003028  0.002390  0.073970
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.005072   0.016647   0.305    0.769
# us          0.935543   0.107872   8.673 5.42e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.03421 on 7 degrees of freedom
# Multiple R-squared:  0.9149,	Adjusted R-squared:  0.9027
# F-statistic: 75.22 on 1 and 7 DF,  p-value: 5.424e-05

p2 <- ggplot(GABA, mapping = aes(x = edlein, y = us)) + geom_point(color = "grey30") + geom_smooth(method = "lm", se = T) + coord_fixed(xlim = c(0, max(GABA[c("edlein", "us")])),
                                                                                                                                        ylim = c(0, max(GABA[c("edlein", "us")])),
                                                                                                                                        expand = F) + cowplot::theme_cowplot() + labs(title = "subclass", subtitle = "Adjusted R-squared:  0.9027\n P = 5.42e-05")

# patchwork::wrap_plots(p1,p2)

# IT
IT2 <- wholeMeta %>% filter(str_detect(cross_area_subclass, "IT")) %>%  group_by(donor, batch, cross_area_cluster) %>% summarise(count = n()) %>% group_by(donor, batch) %>% mutate(per = count /
                                                                                                                                                                                      sum(count)) %>% group_by(batch, cross_area_cluster) %>% summarise(aver = median(per)) %>% pivot_wider(names_from = batch, values_from = aver)

md <- lm(edlein ~ us, IT2)

summary(md)
# Call:
#   lm(formula = edlein ~ us, data = IT2)
#
# Residuals:
#   Min        1Q    Median        3Q       Max
# -0.015214 -0.013542 -0.007150  0.001507  0.049522
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.015861   0.004497   3.527  0.00212 **
#   us          0.642027   0.033431  19.204 2.34e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.0198 on 20 degrees of freedom
# Multiple R-squared:  0.9486,	Adjusted R-squared:  0.946
# F-statistic: 368.8 on 1 and 20 DF,  p-value: 2.341e-14


p3 <- ggplot(IT2, mapping = aes(x = edlein, y = us)) + geom_point(color = "grey30") + geom_smooth(method = "lm", se = T) + coord_fixed(xlim = c(0, max(IT2[c("edlein", "us")])),
                                                                                                                                       ylim = c(0, max(IT2[c("edlein", "us")])),
                                                                                                                                       expand = F) + cowplot::theme_cowplot() + labs(title = "cluster", subtitle = "Adjusted R-squared:  0.946\n P = 2.34e-14")


# GABA
GABA2 <- wholeMeta %>% filter(str_detect(class, "inhibitory")) %>%  group_by(donor, batch, cross_area_cluster) %>% summarise(count = n()) %>% group_by(donor, batch) %>% mutate(per = count /
                                                                                                                                                                                  sum(count)) %>% group_by(batch, cross_area_cluster) %>% summarise(aver = median(per)) %>% pivot_wider(names_from = batch, values_from = aver)

md <- lm(edlein ~ us, GABA2)

summary(md)
# Call:
#   lm(formula = edlein ~ us, data = GABA2)
#
# Residuals:
#   Min         1Q     Median         3Q        Max
# -0.0125987 -0.0017442 -0.0005834  0.0010441  0.0182402
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.0023951  0.0004184   5.725 1.26e-07 ***
#   us          0.7766766  0.0180234  43.093  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.003643 on 93 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.9523,	Adjusted R-squared:  0.9518
# F-statistic:  1857 on 1 and 93 DF,  p-value: < 2.2e-16

p4 <- ggplot(GABA2, mapping = aes(x = edlein, y = us)) + geom_point(color = "grey30") + geom_smooth(method = "lm", se = T) + coord_fixed(xlim = c(0, max(GABA2[c("edlein", "us")])),
                                                                                                                                         ylim = c(0, max(GABA2[c("edlein", "us")])),
                                                                                                                                         expand = F) + cowplot::theme_cowplot() + labs(title = "cluster",subtitle = "Adjusted R-squared:  0.9518\n P = 2.2e-16") + theme(legend.position = "none")

global <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
ggsave("figs1",height = 10,width = 10)

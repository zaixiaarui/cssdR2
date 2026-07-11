rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/input"
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/output"
library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)
set.seed(123)

sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
)

nor_cell_sub_raw <- read_csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

combined_db <- read_csv(
  file.path(input, "sarg/ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

colnames(combined_db) <- c(
  "gene", "type", "subtype", "HMM.category",
  "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank"
)

# -----------------------------
# 2. 样本列名统一管理
# -----------------------------
# 优先使用 sample.csv 中 sample 列与丰度表列名的交集。
# 如果识别不到，则使用手动指定的样本名。
sample_cols <- intersect(sam$sample, names(nor_cell_sub_raw))

if (length(sample_cols) == 0) {
  sample_cols <- names(nor_cell_sub_raw)[names(nor_cell_sub_raw) %in% c(
    "BJ", "CC1", "CC2", "CD1", "CD2", "CQ", "CS", "DBC", "FZ", "HF",
    "HHB1", "HHB2", "JN", "KF", "LZ", "NB", "NJ", "NN1", "NN2", "QD",
    "SSJ1", "SSJ2", "SZ", "WF", "WH", "XA", "XM", "YC", "YX", "ZH"
  )]
}

n_sample <- length(sample_cols)

if (n_sample == 0) {
  stop("没有识别到样本列，请检查 sample.csv 的 sample 列和 normalized_cell.subtype.csv 的列名。")
}

message("Number of sample columns: ", n_sample)
message("Sample columns: ", paste(sample_cols, collapse = ", "))

# 仅保留 subtype 和样本丰度列，避免原始表中已有 type 等注释列导致 join 后出现 type.x/type.y。
nor_cell_sub_abun <- nor_cell_sub_raw %>%
  select(subtype, all_of(sample_cols)) %>%
  group_by(subtype) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# -----------------------------
# 3. 合并 ARG 注释信息
# -----------------------------
# 如果同一个 subtype 在 ARGRANKER_DB 中对应多个 gene，
# 这里保留每个 subtype 的第一条注释，避免 join 后重复扩增丰度表。
combined_db_subtype <- combined_db %>%
  select(-gene) %>%
  distinct(subtype, .keep_all = TRUE)

nor_cell_sub_anno <- nor_cell_sub_abun %>%
  left_join(combined_db_subtype, by = "subtype") %>%
  mutate(
    type = replace_na(type, "others"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  arrange(Rank, subtype)

#4.因子表
factors <- read_csv(
  file.path(input, "factors0527.csv"),
  show_col_types = FALSE
)
head(factors)
head(nor_cell_sub_anno)









#统计分析


按照ARGs总丰度进行分析

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(mlr)
library(randomForest)
library(e1071)

outp <- "outp/factor_ARGs_total_ML_0603"
dir.create(outp, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 整理环境因子 factors
# ============================================================

factor_df <- factors %>%
  mutate(across(everything(), ~na_if(as.character(.x), "/"))) %>%
  mutate(across(
    c(As, Hg, P, Cd, Cr, Pb, N, OM,
      `Annual average temperature`,
      `Annual precipitation`,
      `Annual sunshine hours`,
      `Green area`,
      `Per capita regional GDP`,
      `Total population`,
      longitude, latitude),
    as.numeric
  )) %>%
  mutate(
    ktype = factor(ktype),
    type1 = factor(type1),
    source = factor(source),
    `climate type` = factor(`climate type`)
  ) %>%
  column_to_rownames("sample")

# ============================================================
# 2. 整理 ARG subtype 丰度矩阵
#    行 = sample，列 = ARG subtype
# ============================================================

arg_mat <- nor_cell_sub_anno %>%
  select(subtype, where(is.numeric)) %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat <- arg_mat[, colSums(arg_mat, na.rm = TRUE) > 0]

common_samples <- intersect(rownames(arg_mat), rownames(factor_df))

arg_mat <- arg_mat[common_samples, ]
factor_df <- factor_df[common_samples, ]

# ============================================================
# 3. 提取数值型因子
# ============================================================

env_num <- factor_df %>%
  select(As, Hg, P, Cd, Cr, Pb, N, OM,
         `Annual average temperature`,
         `Annual precipitation`,
         `Annual sunshine hours`,
         `Green area`,
         `Per capita regional GDP`,
         `Total population`,
         longitude, latitude)

# ============================================================
# 4. 计算总 ARG 丰度
# ============================================================

arg_total <- data.frame(
  sample = rownames(arg_mat),
  total_ARG_abundance = rowSums(arg_mat, na.rm = TRUE)
)

min_pos <- min(arg_total$total_ARG_abundance[arg_total$total_ARG_abundance > 0], na.rm = TRUE)

arg_total <- arg_total %>%
  mutate(
    log10_total_ARG_abundance = log10(total_ARG_abundance + min_pos / 2)
  ) %>%
  left_join(
    factor_df %>% rownames_to_column("sample"),
    by = "sample"
  )

write.csv(
  arg_total,
  file.path(outp, "01_total_ARG_abundance_with_factors.csv"),
  row.names = FALSE
)
# ============================================================
# 5. Spearman 相关性：总 ARG 丰度 vs 环境/经济因子
# ============================================================

cor_total <- tibble()

for (i in colnames(env_num)) {
  
  tmp <- data.frame(
    ARG = arg_total$total_ARG_abundance,
    factor_value = env_num[rownames(arg_mat), i]
  ) %>%
    drop_na()
  
  if (nrow(tmp) >= 5 && sd(tmp$ARG) > 0 && sd(tmp$factor_value) > 0) {
    
    test <- cor.test(
      tmp$ARG,
      tmp$factor_value,
      method = "spearman",
      exact = FALSE
    )
    
    cor_total <- bind_rows(
      cor_total,
      data.frame(
        factor = i,
        n = nrow(tmp),
        rho = unname(test$estimate),
        p_value = test$p.value
      )
    )
  }
}

cor_total <- cor_total %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  cor_total,
  file.path(outp, "02_total_ARG_vs_factors_spearman.csv"),
  row.names = FALSE
)

p_cor <- cor_total %>%
  mutate(factor = factor(factor, levels = rev(factor))) %>%
  ggplot(aes(x = factor, y = rho, fill = p_adj < 0.05)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spearman rho",
    fill = "BH-adjusted P < 0.05",
    title = "Correlation between total ARG abundance and factors"
  )
p_cor
ggsave(
  file.path(outp, "02_total_ARG_vs_factors_spearman.pdf"),
  p_cor,
  width = 7,
  height = 5
)

# ============================================================
# 6. 二次项非线性拟合
# ============================================================

quad_summary <- tibble()

for (i in colnames(env_num)) {
  
  tmp <- data.frame(
    ARG = arg_total$log10_total_ARG_abundance,
    factor_value = env_num[rownames(arg_mat), i]
  ) %>%
    drop_na()
  
  if (nrow(tmp) >= 6 && sd(tmp$ARG) > 0 && sd(tmp$factor_value) > 0) {
    
    fit_linear <- lm(ARG ~ factor_value, data = tmp)
    fit_quad <- lm(ARG ~ factor_value + I(factor_value^2), data = tmp)
    
    fit_quad_sum <- summary(fit_quad)
    anova_compare <- anova(fit_linear, fit_quad)
    
    beta1 <- coef(fit_quad)[2]
    beta2 <- coef(fit_quad)[3]
    
    turning_point <- -beta1 / (2 * beta2)
    
    quad_summary <- bind_rows(
      quad_summary,
      data.frame(
        factor = i,
        n = nrow(tmp),
        linear_R2 = summary(fit_linear)$r.squared,
        quadratic_R2 = fit_quad_sum$r.squared,
        quadratic_adj_R2 = fit_quad_sum$adj.r.squared,
        AIC_linear = AIC(fit_linear),
        AIC_quadratic = AIC(fit_quad),
        beta_linear = beta1,
        beta_quadratic = beta2,
        turning_point = turning_point,
        p_linear = coef(fit_quad_sum)[2, 4],
        p_quadratic = coef(fit_quad_sum)[3, 4],
        p_anova_linear_vs_quadratic = anova_compare$`Pr(>F)`[2]
      )
    )
  }
}

quad_summary <- quad_summary %>%
  mutate(
    p_adj_quadratic = p.adjust(p_quadratic, method = "BH"),
    p_adj_anova = p.adjust(p_anova_linear_vs_quadratic, method = "BH"),
    nonlinear_pattern = case_when(
      beta_quadratic < 0 & p_adj_quadratic < 0.05 ~ "inverted_U",
      beta_quadratic > 0 & p_adj_quadratic < 0.05 ~ "U_shape",
      TRUE ~ "not_significant"
    )
  ) %>%
  arrange(p_adj_quadratic)

write.csv(
  quad_summary,
  file.path(outp, "03_total_ARG_nonlinear_quadratic_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 7. 输出每个因子的二次项 + LOESS 拟合图
# ============================================================

pdf(
  file.path(outp, "03_total_ARG_nonlinear_quadratic_LOESS_plots.pdf"),
  width = 6,
  height = 5
)

for (i in colnames(env_num)) {
  
  tmp <- data.frame(
    sample = rownames(arg_mat),
    ARG = arg_total$log10_total_ARG_abundance,
    factor_value = env_num[rownames(arg_mat), i],
    ktype = factor_df[rownames(arg_mat), "ktype"]
  ) %>%
    drop_na()
  
  if (nrow(tmp) >= 6 && sd(tmp$ARG) > 0 && sd(tmp$factor_value) > 0) {
    
    p <- ggplot(tmp, aes(x = factor_value, y = ARG)) +
      geom_point(aes(color = ktype), size = 3) +
      geom_smooth(
        method = "lm",
        formula = y ~ poly(x, 2, raw = TRUE),
        se = TRUE,
        linewidth = 0.8
      ) +
      geom_smooth(
        method = "loess",
        se = FALSE,
        linetype = 2,
        linewidth = 0.8
      ) +
      theme_bw() +
      labs(
        x = i,
        y = "log10(total ARG abundance)",
        color = "ARG group",
        title = paste0("Nonlinear fitting: ", i)
      )
    
    print(p)
  }
}

dev.off()

# ============================================================
# 8. 准备 mlr 数据
# ============================================================

ml_factor_raw <- env_num

# 去掉缺失太多的因子，要求至少 70% 样本有值
factor_keep <- names(which(colSums(!is.na(ml_factor_raw)) >= 0.7 * nrow(ml_factor_raw)))

ml_factor <- ml_factor_raw[, factor_keep, drop = FALSE]

# 去掉零方差因子
factor_sd <- apply(ml_factor, 2, sd, na.rm = TRUE)
factor_keep <- names(factor_sd[factor_sd > 0])

ml_factor <- ml_factor[, factor_keep, drop = FALSE]

# 用 tibble 防止列名被自动修改
ml_data <- bind_cols(
  tibble(
    sample = rownames(arg_mat),
    log10_total_ARG_abundance = arg_total$log10_total_ARG_abundance
  ),
  as_tibble(ml_factor, .name_repair = "minimal")
)

# 中位数填补缺失值
for (i in factor_keep) {
  med_value <- median(ml_data[[i]], na.rm = TRUE)
  ml_data[[i]][is.na(ml_data[[i]])] <- med_value
}

# mlr 不适合使用带空格的变量名，生成安全变量名
factor_name_map <- data.frame(
  original_factor = factor_keep,
  mlr_factor = make.names(factor_keep, unique = TRUE),
  stringsAsFactors = FALSE
)

# 替换环境因子列名
colnames(ml_data)[match(factor_name_map$original_factor, colnames(ml_data))] <- factor_name_map$mlr_factor

# 检查是否还有 NA
match(factor_name_map$original_factor, c("sample", "log10_total_ARG_abundance", factor_name_map$original_factor))

# 去掉 sample，作为 mlr 输入
ml_data_mlr <- ml_data %>%
  select(-sample)

write.csv(
  factor_name_map,
  file.path(outp, "04_mlr_factor_name_map.csv"),
  row.names = FALSE
)

write.csv(
  ml_data,
  file.path(outp, "04_mlr_input_data.csv"),
  row.names = FALSE
)
# ============================================================
# 9. 构建 mlr 回归任务
# ============================================================

task_arg <- makeRegrTask(
  id = "total_ARG_abundance",
  data = ml_data_mlr,
  target = "log10_total_ARG_abundance"
)

# 随机森林
rf_lrn <- makeLearner(
  "regr.randomForest",
  ntree = 1000,
  importance = TRUE
)

# 支持向量机
svm_lrn <- makeLearner(
  "regr.svm",
  kernel = "radial",
  cost = 1,
  gamma = 1 / length(factor_name_map$mlr_factor),
  scale = TRUE
)

# ============================================================
# 10. 交叉验证评估模型性能
# ============================================================

cv_iters <- min(5, nrow(ml_data_mlr))

rdesc <- makeResampleDesc(
  "CV",
  iters = cv_iters
)

set.seed(123)

rf_cv <- resample(
  learner = rf_lrn,
  task = task_arg,
  resampling = rdesc,
  measures = list(rmse, mae, rsq),
  show.info = FALSE
)

set.seed(123)

svm_cv <- resample(
  learner = svm_lrn,
  task = task_arg,
  resampling = rdesc,
  measures = list(rmse, mae, rsq),
  show.info = FALSE
)

model_performance <- data.frame(
  model = c("Random forest", "Support vector machine"),
  RMSE = c(rf_cv$aggr["rmse.test.mean"], svm_cv$aggr["rmse.test.mean"]),
  MAE = c(rf_cv$aggr["mae.test.mean"], svm_cv$aggr["mae.test.mean"]),
  R2 = c(rf_cv$aggr["rsq.test.mean"], svm_cv$aggr["rsq.test.mean"])
)

write.csv(
  model_performance,
  file.path(outp, "05_mlr_RF_SVM_cross_validation_performance.csv"),
  row.names = FALSE
)

# ============================================================
# 11. 训练最终 RF 和 SVM 模型
# ============================================================

set.seed(123)
rf_model <- train(rf_lrn, task_arg)

set.seed(123)
svm_model <- train(svm_lrn, task_arg)

rf_pred <- predict(rf_model, task = task_arg)
svm_pred <- predict(svm_model, task = task_arg)

train_performance <- data.frame(
  model = c("Random forest", "Support vector machine"),
  RMSE = c(
    performance(rf_pred, measures = rmse),
    performance(svm_pred, measures = rmse)
  ),
  MAE = c(
    performance(rf_pred, measures = mae),
    performance(svm_pred, measures = mae)
  ),
  R2 = c(
    performance(rf_pred, measures = rsq),
    performance(svm_pred, measures = rsq)
  )
)

write.csv(
  train_performance,
  file.path(outp, "06_mlr_RF_SVM_training_performance.csv"),
  row.names = FALSE
)

# ============================================================
# 12. RF 内部变量重要性
# ============================================================

rf_inner <- getLearnerModel(rf_model)

rf_internal_importance <- randomForest::importance(rf_inner) %>%
  as.data.frame() %>%
  rownames_to_column("mlr_factor") %>%
  left_join(factor_name_map, by = "mlr_factor") %>%
  arrange(desc(`%IncMSE`))

write.csv(
  rf_internal_importance,
  file.path(outp, "07_RF_internal_variable_importance.csv"),
  row.names = FALSE
)

p_rf_internal <- rf_internal_importance %>%
  slice_head(n = 20) %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "%IncMSE",
    title = "Random forest internal variable importance"
  )

ggsave(
  file.path(outp, "07_RF_internal_variable_importance.pdf"),
  p_rf_internal,
  width = 7,
  height = 5
)

# ============================================================
# 13. RF permutation importance
# ============================================================

set.seed(123)

rf_base_pred <- predict(rf_model, task = task_arg)
rf_base_rmse <- performance(rf_base_pred, measures = rmse)

rf_perm_importance <- tibble()

for (j in factor_name_map$mlr_factor) {
  
  perm_values <- c()
  
  for (k in 1:100) {
    
    tmp_data <- ml_data_mlr
    tmp_data[[j]] <- sample(tmp_data[[j]])
    
    tmp_task <- makeRegrTask(
      id = "tmp_rf_perm",
      data = tmp_data,
      target = "log10_total_ARG_abundance"
    )
    
    tmp_pred <- predict(rf_model, task = tmp_task)
    tmp_rmse <- performance(tmp_pred, measures = rmse)
    
    perm_values <- c(perm_values, tmp_rmse - rf_base_rmse)
  }
  
  rf_perm_importance <- bind_rows(
    rf_perm_importance,
    data.frame(
      mlr_factor = j,
      importance_mean = mean(perm_values),
      importance_sd = sd(perm_values),
      importance_positive_rate = mean(perm_values > 0)
    )
  )
}

rf_perm_importance <- rf_perm_importance %>%
  left_join(factor_name_map, by = "mlr_factor") %>%
  arrange(desc(importance_mean))

write.csv(
  rf_perm_importance,
  file.path(outp, "08_RF_permutation_importance.csv"),
  row.names = FALSE
)

p_rf_perm <- rf_perm_importance %>%
  slice_head(n = 20) %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = importance_mean)) +
  geom_col() +
  geom_errorbar(
    aes(
      ymin = importance_mean - importance_sd,
      ymax = importance_mean + importance_sd
    ),
    width = 0.2
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Increase in RMSE after permutation",
    title = "Random forest permutation importance"
  )

ggsave(
  file.path(outp, "08_RF_permutation_importance.pdf"),
  p_rf_perm,
  width = 7,
  height = 5
)

# ============================================================
# 14. SVM permutation importance
# ============================================================

set.seed(123)

svm_base_pred <- predict(svm_model, task = task_arg)
svm_base_rmse <- performance(svm_base_pred, measures = rmse)

svm_perm_importance <- tibble()

for (j in factor_name_map$mlr_factor) {
  
  perm_values <- c()
  
  for (k in 1:100) {
    
    tmp_data <- ml_data_mlr
    tmp_data[[j]] <- sample(tmp_data[[j]])
    
    tmp_task <- makeRegrTask(
      id = "tmp_svm_perm",
      data = tmp_data,
      target = "log10_total_ARG_abundance"
    )
    
    tmp_pred <- predict(svm_model, task = tmp_task)
    tmp_rmse <- performance(tmp_pred, measures = rmse)
    
    perm_values <- c(perm_values, tmp_rmse - svm_base_rmse)
  }
  
  svm_perm_importance <- bind_rows(
    svm_perm_importance,
    data.frame(
      mlr_factor = j,
      importance_mean = mean(perm_values),
      importance_sd = sd(perm_values),
      importance_positive_rate = mean(perm_values > 0)
    )
  )
}

svm_perm_importance <- svm_perm_importance %>%
  left_join(factor_name_map, by = "mlr_factor") %>%
  arrange(desc(importance_mean))

write.csv(
  svm_perm_importance,
  file.path(outp, "09_SVM_permutation_importance.csv"),
  row.names = FALSE
)

p_svm_perm <- svm_perm_importance %>%
  slice_head(n = 20) %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = importance_mean)) +
  geom_col() +
  geom_errorbar(
    aes(
      ymin = importance_mean - importance_sd,
      ymax = importance_mean + importance_sd
    ),
    width = 0.2
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Increase in RMSE after permutation",
    title = "SVM permutation importance"
  )

ggsave(
  file.path(outp, "09_SVM_permutation_importance.pdf"),
  p_svm_perm,
  width = 7,
  height = 5
)

# ============================================================
# 15. 综合 RF 和 SVM 重要性排序
# ============================================================

rf_rank <- rf_perm_importance %>%
  select(original_factor, RF_importance = importance_mean) %>%
  mutate(RF_rank = rank(-RF_importance, ties.method = "average"))

svm_rank <- svm_perm_importance %>%
  select(original_factor, SVM_importance = importance_mean) %>%
  mutate(SVM_rank = rank(-SVM_importance, ties.method = "average"))

ml_importance_summary <- rf_rank %>%
  full_join(svm_rank, by = "original_factor") %>%
  mutate(
    mean_rank = rowMeans(select(., RF_rank, SVM_rank), na.rm = TRUE),
    mean_importance = rowMeans(select(., RF_importance, SVM_importance), na.rm = TRUE)
  ) %>%
  arrange(mean_rank)

write.csv(
  ml_importance_summary,
  file.path(outp, "10_RF_SVM_integrated_important_factors.csv"),
  row.names = FALSE
)

p_integrated <- ml_importance_summary %>%
  slice_head(n = 20) %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = mean_rank)) +
  geom_col() +
  coord_flip() +
  scale_y_reverse() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Mean rank of RF and SVM",
    title = "Integrated important factors selected by RF and SVM"
  )

ggsave(
  file.path(outp, "10_RF_SVM_integrated_important_factors.pdf"),
  p_integrated,
  width = 7,
  height = 5
)




















按照所有ARGs进行分析

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(mlr)
library(randomForest)
library(e1071)

outp <- "outp/factor_ARGs_profile"
dir.create(outp, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 整理环境因子
# ============================================================

factor_df <- factors %>%
  mutate(across(everything(), ~na_if(as.character(.x), "/"))) %>%
  mutate(across(
    c(As, Hg, P, Cd, Cr, Pb, N, OM,
      `Annual average temperature`,
      `Annual precipitation`,
      `Annual sunshine hours`,
      `Green area`,
      `Per capita regional GDP`,
      `Total population`,
      longitude, latitude),
    as.numeric
  )) %>%
  mutate(
    ktype = factor(ktype),
    type1 = factor(type1),
    source = factor(source),
    `climate type` = factor(`climate type`)
  ) %>%
  column_to_rownames("sample")

# ============================================================
# 2. 整理整体 ARG subtype 丰度矩阵
# ============================================================

arg_mat <- nor_cell_sub_anno %>%
  select(subtype, where(is.numeric)) %>%
  group_by(subtype) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat <- arg_mat[, colSums(arg_mat, na.rm = TRUE) > 0]

common_samples <- intersect(rownames(arg_mat), rownames(factor_df))

arg_mat <- arg_mat[common_samples, ]
factor_df <- factor_df[common_samples, ]

# ============================================================
# 3. 提取数值型因子
# ============================================================

env_num <- factor_df %>%
  select(As, Hg, P, Cd, Cr, Pb, N, OM,
         `Annual average temperature`,
         `Annual precipitation`,
         `Annual sunshine hours`,
         `Green area`,
         `Per capita regional GDP`,
         `Total population`,
         longitude, latitude)

# 去掉缺失太多的因子
env_num <- env_num[, colSums(!is.na(env_num)) >= 0.7 * nrow(env_num), drop = FALSE]

# 去掉零方差因子
env_sd <- apply(env_num, 2, sd, na.rm = TRUE)
env_num <- env_num[, env_sd > 0, drop = FALSE]

# 中位数填补缺失
for (i in colnames(env_num)) {
  env_num[[i]][is.na(env_num[[i]])] <- median(env_num[[i]], na.rm = TRUE)
}

write.csv(
  env_num %>% rownames_to_column("sample"),
  file.path(outp, "01_factor_matrix_used.csv"),
  row.names = FALSE
)

write.csv(
  arg_mat %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_subtype_matrix_used.csv"),
  row.names = FALSE
)


# ============================================================
# 4. Bray-Curtis 距离
# ============================================================

arg_bray <- vegdist(arg_mat, method = "bray")

# ============================================================
# 5. 单因子 PERMANOVA
# ============================================================

permanova_single <- tibble()

for (i in colnames(env_num)) {
  
  tmp_env <- data.frame(
    factor_value = env_num[[i]]
  )
  
  rownames(tmp_env) <- rownames(env_num)
  
  test <- adonis2(
    arg_bray ~ factor_value,
    data = tmp_env,
    permutations = 999
  )
  
  permanova_single <- bind_rows(
    permanova_single,
    data.frame(
      factor = i,
      R2 = test$R2[1],
      F_value = test$F[1],
      p_value = test$`Pr(>F)`[1]
    )
  )
}

permanova_single <- permanova_single %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  permanova_single,
  file.path(outp, "02_ARG_profile_single_factor_PERMANOVA.csv"),
  row.names = FALSE
)

p_permanova <- permanova_single %>%
  mutate(factor = factor(factor, levels = rev(factor))) %>%
  ggplot(aes(x = factor, y = R2, fill = p_adj < 0.05)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "PERMANOVA R2",
    fill = "BH-adjusted P < 0.05",
    title = "Single-factor effects on overall ARG composition"
  )

ggsave(
  file.path(outp, "02_ARG_profile_single_factor_PERMANOVA.pdf"),
  p_permanova,
  width = 7,
  height = 5
)


# ============================================================
# 6. db-RDA 分析
# ============================================================

env_scale <- as.data.frame(scale(env_num))

dbrda_full <- capscale(arg_bray ~ ., data = env_scale)

write.csv(
  as.data.frame(anova.cca(dbrda_full, permutations = 999)),
  file.path(outp, "03_dbRDA_overall_test.csv")
)

write.csv(
  as.data.frame(anova.cca(dbrda_full, by = "axis", permutations = 999)),
  file.path(outp, "03_dbRDA_axis_test.csv")
)

write.csv(
  as.data.frame(anova.cca(dbrda_full, by = "term", permutations = 999)),
  file.path(outp, "03_dbRDA_factor_test.csv")
)

write.csv(
  data.frame(
    factor = names(vif.cca(dbrda_full)),
    VIF = as.numeric(vif.cca(dbrda_full))
  ),
  file.path(outp, "03_dbRDA_VIF.csv"),
  row.names = FALSE
)


# ============================================================
# 7. db-RDA 作图
# ============================================================

site_score <- scores(dbrda_full, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(factor_df %>% rownames_to_column("sample"), by = "sample")

env_score <- scores(dbrda_full, display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column("factor")

dbrda_imp <- summary(dbrda_full)$cont$importance

dbrda_exp1 <- ifelse(ncol(dbrda_imp) >= 1, dbrda_imp[2, 1] * 100, 0)
dbrda_exp2 <- ifelse(ncol(dbrda_imp) >= 2, dbrda_imp[2, 2] * 100, 0)

p_dbrda <- ggplot(site_score, aes(CAP1, CAP2, color = ktype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = ktype), linetype = 2) +
  geom_segment(
    data = env_score,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    arrow = arrow(length = unit(0.25, "cm")),
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = env_score,
    aes(CAP1, CAP2, label = factor),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw() +
  labs(
    x = paste0("CAP1 (", round(dbrda_exp1, 2), "%)"),
    y = paste0("CAP2 (", round(dbrda_exp2, 2), "%)"),
    color = "ARG group",
    title = "db-RDA of overall ARG composition"
  )

ggsave(
  file.path(outp, "03_dbRDA_ARG_profile_factors.pdf"),
  p_dbrda,
  width = 7,
  height = 6
)


# ============================================================
# 8. PCoA + envfit
# ============================================================

pcoa_arg <- cmdscale(arg_bray, eig = TRUE, k = 2)

pcoa_site <- data.frame(
  sample = rownames(arg_mat),
  PCoA1 = pcoa_arg$points[, 1],
  PCoA2 = pcoa_arg$points[, 2]
) %>%
  left_join(factor_df %>% rownames_to_column("sample"), by = "sample")

envfit_res <- envfit(
  pcoa_arg$points,
  env_scale[rownames(pcoa_arg$points), ],
  permutations = 999
)

envfit_df <- as.data.frame(scores(envfit_res, display = "vectors")) %>%
  rownames_to_column("factor") %>%
  mutate(
    r2 = envfit_res$vectors$r,
    p_value = envfit_res$vectors$pvals,
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  envfit_df,
  file.path(outp, "04_envfit_ARG_profile_factors.csv"),
  row.names = FALSE
)

p_envfit <- envfit_df %>%
  mutate(factor = factor(factor, levels = rev(factor))) %>%
  ggplot(aes(x = factor, y = r2, fill = p_adj < 0.05)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "envfit r2",
    fill = "BH-adjusted P < 0.05",
    title = "Environmental fitting for overall ARG composition"
  )

ggsave(
  file.path(outp, "04_envfit_ARG_profile_factors.pdf"),
  p_envfit,
  width = 7,
  height = 5
)


# ============================================================
# 9. PCoA 图
# ============================================================

eig_value <- pcoa_arg$eig
eig_percent <- eig_value / sum(eig_value[eig_value > 0]) * 100

p_pcoa <- ggplot(pcoa_site, aes(PCoA1, PCoA2, color = ktype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = ktype), linetype = 2) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", round(eig_percent[1], 2), "%)"),
    y = paste0("PCoA2 (", round(eig_percent[2], 2), "%)"),
    color = "ARG group",
    title = "PCoA of overall ARG composition"
  )

ggsave(
  file.path(outp, "05_PCoA_ARG_profile_ktype.pdf"),
  p_pcoa,
  width = 6,
  height = 5
)


# ============================================================
# 10. 准备 mlr 数据：预测 ARGs 组成主轴
# ============================================================

ml_factor <- env_num

factor_keep <- colnames(ml_factor)

ml_data <- bind_cols(
  tibble(
    sample = rownames(arg_mat),
    PCoA1 = pcoa_arg$points[, 1],
    PCoA2 = pcoa_arg$points[, 2]
  ),
  as_tibble(ml_factor, .name_repair = "minimal")
)

factor_name_map <- data.frame(
  original_factor = factor_keep,
  mlr_factor = make.names(factor_keep, unique = TRUE),
  stringsAsFactors = FALSE
)

colnames(ml_data)[match(factor_name_map$original_factor, colnames(ml_data))] <- factor_name_map$mlr_factor

ml_data_pcoa1 <- ml_data %>%
  select(-sample, -PCoA2)

ml_data_pcoa2 <- ml_data %>%
  select(-sample, -PCoA1)

write.csv(
  factor_name_map,
  file.path(outp, "06_mlr_factor_name_map.csv"),
  row.names = FALSE
)

write.csv(
  ml_data,
  file.path(outp, "06_mlr_input_ARG_profile_PCoA.csv"),
  row.names = FALSE
)

# ============================================================
# 11. PCoA1 任务
# ============================================================

task_pcoa1 <- makeRegrTask(
  id = "ARG_profile_PCoA1",
  data = ml_data_pcoa1,
  target = "PCoA1"
)

rf_lrn <- makeLearner(
  "regr.randomForest",
  ntree = 1000,
  importance = TRUE
)

svm_lrn <- makeLearner(
  "regr.svm",
  kernel = "radial",
  cost = 1,
  gamma = 1 / length(factor_name_map$mlr_factor),
  scale = TRUE
)

set.seed(123)
rf_model_pcoa1 <- train(rf_lrn, task_pcoa1)

set.seed(123)
svm_model_pcoa1 <- train(svm_lrn, task_pcoa1)

# ============================================================
# 12. PCoA2 任务
# ============================================================

task_pcoa2 <- makeRegrTask(
  id = "ARG_profile_PCoA2",
  data = ml_data_pcoa2,
  target = "PCoA2"
)

set.seed(123)
rf_model_pcoa2 <- train(rf_lrn, task_pcoa2)

set.seed(123)
svm_model_pcoa2 <- train(svm_lrn, task_pcoa2)


# ============================================================
# 13. RF permutation importance for PCoA1
# ============================================================

rf_base_pred_pcoa1 <- predict(rf_model_pcoa1, task = task_pcoa1)
rf_base_rmse_pcoa1 <- performance(rf_base_pred_pcoa1, measures = rmse)

rf_perm_pcoa1 <- tibble()

for (j in factor_name_map$mlr_factor) {
  
  perm_values <- c()
  
  for (k in 1:100) {
    
    tmp_data <- ml_data_pcoa1
    tmp_data[[j]] <- sample(tmp_data[[j]])
    
    tmp_task <- makeRegrTask(
      id = "tmp_rf_pcoa1",
      data = tmp_data,
      target = "PCoA1"
    )
    
    tmp_pred <- predict(rf_model_pcoa1, task = tmp_task)
    tmp_rmse <- performance(tmp_pred, measures = rmse)
    
    perm_values <- c(perm_values, tmp_rmse - rf_base_rmse_pcoa1)
  }
  
  rf_perm_pcoa1 <- bind_rows(
    rf_perm_pcoa1,
    data.frame(
      mlr_factor = j,
      RF_PCoA1_importance = mean(perm_values)
    )
  )
}

# ============================================================
# 14. SVM permutation importance for PCoA1
# ============================================================

svm_base_pred_pcoa1 <- predict(svm_model_pcoa1, task = task_pcoa1)
svm_base_rmse_pcoa1 <- performance(svm_base_pred_pcoa1, measures = rmse)

svm_perm_pcoa1 <- tibble()

for (j in factor_name_map$mlr_factor) {
  
  perm_values <- c()
  
  for (k in 1:100) {
    
    tmp_data <- ml_data_pcoa1
    tmp_data[[j]] <- sample(tmp_data[[j]])
    
    tmp_task <- makeRegrTask(
      id = "tmp_svm_pcoa1",
      data = tmp_data,
      target = "PCoA1"
    )
    
    tmp_pred <- predict(svm_model_pcoa1, task = tmp_task)
    tmp_rmse <- performance(tmp_pred, measures = rmse)
    
    perm_values <- c(perm_values, tmp_rmse - svm_base_rmse_pcoa1)
  }
  
  svm_perm_pcoa1 <- bind_rows(
    svm_perm_pcoa1,
    data.frame(
      mlr_factor = j,
      SVM_PCoA1_importance = mean(perm_values)
    )
  )
}

# ============================================================
# 15. RF permutation importance for PCoA2
# ============================================================

rf_base_pred_pcoa2 <- predict(rf_model_pcoa2, task = task_pcoa2)
rf_base_rmse_pcoa2 <- performance(rf_base_pred_pcoa2, measures = rmse)

rf_perm_pcoa2 <- tibble()

for (j in factor_name_map$mlr_factor) {
  
  perm_values <- c()
  
  for (k in 1:100) {
    
    tmp_data <- ml_data_pcoa2
    tmp_data[[j]] <- sample(tmp_data[[j]])
    
    tmp_task <- makeRegrTask(
      id = "tmp_rf_pcoa2",
      data = tmp_data,
      target = "PCoA2"
    )
    
    tmp_pred <- predict(rf_model_pcoa2, task = tmp_task)
    tmp_rmse <- performance(tmp_pred, measures = rmse)
    
    perm_values <- c(perm_values, tmp_rmse - rf_base_rmse_pcoa2)
  }
  
  rf_perm_pcoa2 <- bind_rows(
    rf_perm_pcoa2,
    data.frame(
      mlr_factor = j,
      RF_PCoA2_importance = mean(perm_values)
    )
  )
}

# ============================================================
# 16. SVM permutation importance for PCoA2
# ============================================================

svm_base_pred_pcoa2 <- predict(svm_model_pcoa2, task = task_pcoa2)
svm_base_rmse_pcoa2 <- performance(svm_base_pred_pcoa2, measures = rmse)

svm_perm_pcoa2 <- tibble()

for (j in factor_name_map$mlr_factor) {
  
  perm_values <- c()
  
  for (k in 1:100) {
    
    tmp_data <- ml_data_pcoa2
    tmp_data[[j]] <- sample(tmp_data[[j]])
    
    tmp_task <- makeRegrTask(
      id = "tmp_svm_pcoa2",
      data = tmp_data,
      target = "PCoA2"
    )
    
    tmp_pred <- predict(svm_model_pcoa2, task = tmp_task)
    tmp_rmse <- performance(tmp_pred, measures = rmse)
    
    perm_values <- c(perm_values, tmp_rmse - svm_base_rmse_pcoa2)
  }
  
  svm_perm_pcoa2 <- bind_rows(
    svm_perm_pcoa2,
    data.frame(
      mlr_factor = j,
      SVM_PCoA2_importance = mean(perm_values)
    )
  )
}


# ============================================================
# 17. 综合 PCoA1 和 PCoA2 的变量重要性
# ============================================================

ml_importance_ARG_profile <- rf_perm_pcoa1 %>%
  full_join(svm_perm_pcoa1, by = "mlr_factor") %>%
  full_join(rf_perm_pcoa2, by = "mlr_factor") %>%
  full_join(svm_perm_pcoa2, by = "mlr_factor") %>%
  left_join(factor_name_map, by = "mlr_factor") %>%
  mutate(
    mean_importance = rowMeans(
      select(
        .,
        RF_PCoA1_importance,
        SVM_PCoA1_importance,
        RF_PCoA2_importance,
        SVM_PCoA2_importance
      ),
      na.rm = TRUE
    ),
    rank_importance = rank(-mean_importance, ties.method = "average")
  ) %>%
  arrange(rank_importance)

write.csv(
  ml_importance_ARG_profile,
  file.path(outp, "07_RF_SVM_integrated_importance_ARG_profile.csv"),
  row.names = FALSE
)

p_ml_profile <- ml_importance_ARG_profile %>%
  slice_head(n = 20) %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = mean_importance)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Mean increase in RMSE after permutation",
    title = "RF and SVM importance for overall ARG composition"
  )

ggsave(
  file.path(outp, "07_RF_SVM_integrated_importance_ARG_profile.pdf"),
  p_ml_profile,
  width = 7,
  height = 5
)













按照环境、地理和经济分类因子

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(mlr)
library(randomForest)
library(e1071)

outp <- "outp/factor_ARGs_profile_grouped"
dir.create(outp, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 因子分类
# ============================================================

environmental_vars <- c("As", "Hg", "P", "Cd", "Cr", "Pb", "N", "OM")

geo_climate_num_vars <- c(
  "Annual_average_temperature",
  "Annual_precipitation",
  "Annual_sunshine_hours",
  "longitude",
  "latitude"
)

geo_climate_cat_vars <- c("climate_type")

economic_vars <- c(
  "Green_area",
  "Per_capita_regional_GDP",
  "Total_population"
)

factor_group_map_all <- tibble(
  factor = c(
    environmental_vars,
    geo_climate_num_vars,
    geo_climate_cat_vars,
    economic_vars
  ),
  original_factor = c(
    environmental_vars,
    "Annual average temperature",
    "Annual precipitation",
    "Annual sunshine hours",
    "longitude",
    "latitude",
    "climate type",
    "Green area",
    "Per capita regional GDP",
    "Total population"
  ),
  group = c(
    rep("Environmental factors", length(environmental_vars)),
    rep("Geo-climatic factors", length(geo_climate_num_vars)),
    rep("Geo-climatic factors", length(geo_climate_cat_vars)),
    rep("Economic factors", length(economic_vars))
  )
)

# ============================================================
# 2. 整理 factors
# ============================================================

factor_df <- factors %>%
  mutate(across(everything(), ~na_if(as.character(.x), "/"))) %>%
  rename(
    climate_type = `climate type`,
    Annual_average_temperature = `Annual average temperature`,
    Annual_precipitation = `Annual precipitation`,
    Annual_sunshine_hours = `Annual sunshine hours`,
    Green_area = `Green area`,
    Per_capita_regional_GDP = `Per capita regional GDP`,
    Total_population = `Total population`
  ) %>%
  mutate(across(
    all_of(c(environmental_vars, geo_climate_num_vars, economic_vars)),
    as.numeric
  )) %>%
  mutate(
    climate_type = replace_na(climate_type, "Unknown"),
    climate_type = factor(climate_type),
    ktype = factor(ktype)
  ) %>%
  column_to_rownames("sample")

# ============================================================
# 3. 整理整体 ARG subtype 丰度矩阵
# ============================================================

sample_cols <- intersect(colnames(nor_cell_sub_anno), rownames(factor_df))

arg_mat <- nor_cell_sub_anno %>%
  select(subtype, all_of(sample_cols)) %>%
  group_by(subtype) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat[is.na(arg_mat)] <- 0
arg_mat <- arg_mat[, colSums(arg_mat, na.rm = TRUE) > 0, drop = FALSE]

common_samples <- intersect(rownames(arg_mat), rownames(factor_df))

arg_mat <- arg_mat[common_samples, , drop = FALSE]
factor_df <- factor_df[common_samples, , drop = FALSE]

# ============================================================
# 4. ARGs 标准化
#    arg_rel 用于 Bray-Curtis
#    arg_hel 用于 varpart / RDA 类分析
# ============================================================

arg_rel <- decostand(arg_mat, method = "total")
arg_rel[is.na(arg_rel)] <- 0

arg_hel <- decostand(arg_rel, method = "hellinger")
arg_hel[is.na(arg_hel)] <- 0

arg_bray <- vegdist(arg_rel, method = "bray")

write.csv(
  arg_rel %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_relative_abundance_matrix.csv"),
  row.names = FALSE
)

# ============================================================
# 5. 整理数值型因子：去除缺失过多和零方差变量
# ============================================================

all_num_vars <- c(environmental_vars, geo_climate_num_vars, economic_vars)

env_num_raw <- factor_df %>%
  select(any_of(all_num_vars))

num_keep <- names(which(colSums(!is.na(env_num_raw)) >= 0.7 * nrow(env_num_raw)))

env_num <- env_num_raw[, num_keep, drop = FALSE]

num_sd <- apply(env_num, 2, sd, na.rm = TRUE)
num_keep <- names(num_sd[num_sd > 0])

env_num <- env_num[, num_keep, drop = FALSE]

for (i in colnames(env_num)) {
  env_num[[i]][is.na(env_num[[i]])] <- median(env_num[[i]], na.rm = TRUE)
}

env_num_scale <- as.data.frame(scale(env_num))
rownames(env_num_scale) <- rownames(env_num)

# climate type 保留为分类变量
geo_cat <- factor_df %>%
  select(any_of(geo_climate_cat_vars))

geo_cat$climate_type <- droplevels(factor(geo_cat$climate_type))

if (n_distinct(geo_cat$climate_type) > 1) {
  cat_keep <- "climate_type"
} else {
  cat_keep <- character(0)
}

model_df <- cbind(
  env_num_scale,
  geo_cat[, cat_keep, drop = FALSE]
)

rownames(model_df) <- rownames(factor_df)

factor_group_map <- factor_group_map_all %>%
  filter(factor %in% colnames(model_df))

write.csv(
  factor_group_map,
  file.path(outp, "01_factor_group_map_used.csv"),
  row.names = FALSE
)

write.csv(
  model_df %>% rownames_to_column("sample"),
  file.path(outp, "01_factor_matrix_scaled_used.csv"),
  row.names = FALSE
)

# ============================================================
# 6. 单因子 PERMANOVA
# ============================================================

permanova_single <- tibble()

for (i in colnames(model_df)) {
  
  tmp_env <- model_df[, i, drop = FALSE]
  colnames(tmp_env) <- "factor_value"
  
  if (is.numeric(tmp_env$factor_value)) {
    valid_factor <- sd(tmp_env$factor_value, na.rm = TRUE) > 0
  } else {
    valid_factor <- n_distinct(tmp_env$factor_value) > 1
  }
  
  if (valid_factor) {
    
    test <- adonis2(
      arg_bray ~ factor_value,
      data = tmp_env,
      permutations = 999
    )
    
    permanova_single <- bind_rows(
      permanova_single,
      data.frame(
        factor = i,
        R2 = test$R2[1],
        F_value = test$F[1],
        p_value = test$`Pr(>F)`[1]
      )
    )
  }
}

permanova_single <- permanova_single %>%
  left_join(factor_group_map, by = "factor") %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  permanova_single,
  file.path(outp, "02_single_factor_PERMANOVA_grouped.csv"),
  row.names = FALSE
)

p_single <- permanova_single %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = R2, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "PERMANOVA R2",
    fill = "Factor group",
    title = "Single-factor effects on overall ARG composition"
  )

ggsave(
  file.path(outp, "02_single_factor_PERMANOVA_grouped.pdf"),
  p_single,
  width = 7,
  height = 5
)

# ============================================================
# 7. 分组 PERMANOVA
# ============================================================

group_permanova <- tibble()

env_data <- model_df[, intersect(environmental_vars, colnames(model_df)), drop = FALSE]

if (ncol(env_data) > 0) {
  test_env <- adonis2(arg_bray ~ ., data = env_data, permutations = 999)
  
  group_permanova <- bind_rows(
    group_permanova,
    data.frame(
      group = "Environmental factors",
      n_factor = ncol(env_data),
      R2 = test_env$R2[1],
      F_value = test_env$F[1],
      p_value = test_env$`Pr(>F)`[1]
    )
  )
}

geo_data <- model_df[, intersect(c(geo_climate_num_vars, geo_climate_cat_vars), colnames(model_df)), drop = FALSE]

if (ncol(geo_data) > 0) {
  test_geo <- adonis2(arg_bray ~ ., data = geo_data, permutations = 999)
  
  group_permanova <- bind_rows(
    group_permanova,
    data.frame(
      group = "Geo-climatic factors",
      n_factor = ncol(geo_data),
      R2 = test_geo$R2[1],
      F_value = test_geo$F[1],
      p_value = test_geo$`Pr(>F)`[1]
    )
  )
}

econ_data <- model_df[, intersect(economic_vars, colnames(model_df)), drop = FALSE]

if (ncol(econ_data) > 0) {
  test_econ <- adonis2(arg_bray ~ ., data = econ_data, permutations = 999)
  
  group_permanova <- bind_rows(
    group_permanova,
    data.frame(
      group = "Economic factors",
      n_factor = ncol(econ_data),
      R2 = test_econ$R2[1],
      F_value = test_econ$F[1],
      p_value = test_econ$`Pr(>F)`[1]
    )
  )
}

group_permanova <- group_permanova %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(desc(R2))

write.csv(
  group_permanova,
  file.path(outp, "03_factor_group_PERMANOVA.csv"),
  row.names = FALSE
)

p_group <- group_permanova %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = group, y = R2, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = significance), hjust = -0.2, size = 5) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "PERMANOVA R2",
    title = "Grouped effects on overall ARG composition"
  )

ggsave(
  file.path(outp, "03_factor_group_PERMANOVA.pdf"),
  p_group,
  width = 6,
  height = 4
)

# ============================================================
# 8. 全部因子共同作用
# ============================================================

permanova_full <- adonis2(
  arg_bray ~ .,
  data = model_df,
  permutations = 999
)

permanova_term <- adonis2(
  arg_bray ~ .,
  data = model_df,
  permutations = 999,
  by = "term"
)

write.csv(
  as.data.frame(permanova_full),
  file.path(outp, "04_full_model_PERMANOVA.csv")
)

write.csv(
  as.data.frame(permanova_term),
  file.path(outp, "04_full_model_PERMANOVA_by_term.csv")
)

# ============================================================
# 9. varpart 方差分解
# ============================================================

climate_dummy <- model.matrix(~ climate_type, data = model_df)

if (ncol(climate_dummy) > 1) {
  climate_dummy <- climate_dummy[, -1, drop = FALSE]
} else {
  climate_dummy <- matrix(nrow = nrow(model_df), ncol = 0)
}

rownames(climate_dummy) <- rownames(model_df)

env_X <- as.data.frame(
  model_df[, intersect(environmental_vars, colnames(model_df)), drop = FALSE]
)

geo_X <- cbind(
  as.data.frame(model_df[, intersect(geo_climate_num_vars, colnames(model_df)), drop = FALSE]),
  as.data.frame(climate_dummy)
)

econ_X <- as.data.frame(
  model_df[, intersect(economic_vars, colnames(model_df)), drop = FALSE]
)

vp <- varpart(
  arg_hel,
  env_X,
  geo_X,
  econ_X
)

capture.output(
  vp,
  file = file.path(outp, "05_variation_partitioning_ARG_profile.txt")
)

pdf(
  file.path(outp, "05_variation_partitioning_ARG_profile.pdf"),
  width = 6,
  height = 5
)

plot(
  vp,
  Xnames = c(
    "Environmental",
    "Geo-climatic",
    "Economic"
  )
)

dev.off()


# ============================================================
# 10. db-RDA
# ============================================================

dbrda_full <- capscale(
  arg_bray ~ .,
  data = model_df
)

write.csv(
  as.data.frame(anova.cca(dbrda_full, permutations = 999)),
  file.path(outp, "06_dbRDA_overall_test.csv")
)

write.csv(
  as.data.frame(anova.cca(dbrda_full, by = "term", permutations = 999)),
  file.path(outp, "06_dbRDA_factor_test.csv")
)

site_score <- scores(dbrda_full, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(factor_df %>% rownames_to_column("sample"), by = "sample")

env_score <- scores(dbrda_full, display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column("factor") %>%
  left_join(factor_group_map, by = "factor")

dbrda_imp <- summary(dbrda_full)$cont$importance

dbrda_exp1 <- ifelse(ncol(dbrda_imp) >= 1, dbrda_imp[2, 1] * 100, 0)
dbrda_exp2 <- ifelse(ncol(dbrda_imp) >= 2, dbrda_imp[2, 2] * 100, 0)

p_dbrda <- ggplot(site_score, aes(CAP1, CAP2, color = ktype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = ktype), linetype = 2) +
  geom_segment(
    data = env_score,
    aes(
      x = 0,
      y = 0,
      xend = CAP1,
      yend = CAP2,
      linetype = group
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = env_score,
    aes(CAP1, CAP2, label = original_factor),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw() +
  labs(
    x = paste0("CAP1 (", round(dbrda_exp1, 2), "%)"),
    y = paste0("CAP2 (", round(dbrda_exp2, 2), "%)"),
    color = "ARG group",
    linetype = "Factor group",
    title = "db-RDA of overall ARG composition by grouped factors"
  )

ggsave(
  file.path(outp, "06_dbRDA_ARG_profile_grouped_factors.pdf"),
  p_dbrda,
  width = 7,
  height = 6
)

# ============================================================
# 11. PCoA + envfit
# ============================================================

pcoa_arg <- cmdscale(arg_bray, eig = TRUE, k = 2)

pcoa_site <- data.frame(
  sample = rownames(arg_mat),
  PCoA1 = pcoa_arg$points[, 1],
  PCoA2 = pcoa_arg$points[, 2]
) %>%
  left_join(factor_df %>% rownames_to_column("sample"), by = "sample")

envfit_res <- envfit(
  pcoa_arg$points,
  model_df,
  permutations = 999
)

envfit_vec <- as.data.frame(scores(envfit_res, display = "vectors")) %>%
  rownames_to_column("factor") %>%
  mutate(
    r2 = envfit_res$vectors$r,
    p_value = envfit_res$vectors$pvals,
    variable_type = "numeric"
  )

if (!is.null(envfit_res$factors)) {
  
  envfit_fac <- data.frame(
    factor = names(envfit_res$factors$r),
    r2 = as.numeric(envfit_res$factors$r),
    p_value = as.numeric(envfit_res$factors$pvals),
    variable_type = "categorical"
  )
  
} else {
  
  envfit_fac <- tibble(
    factor = character(),
    r2 = numeric(),
    p_value = numeric(),
    variable_type = character()
  )
}

envfit_df <- bind_rows(envfit_vec, envfit_fac) %>%
  left_join(factor_group_map, by = "factor") %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  envfit_df,
  file.path(outp, "07_envfit_ARG_profile_grouped_factors.csv"),
  row.names = FALSE
)

p_envfit <- envfit_df %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = r2, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "envfit r2",
    fill = "Factor group",
    title = "Environmental fitting for overall ARG composition"
  )

ggsave(
  file.path(outp, "07_envfit_ARG_profile_grouped_factors.pdf"),
  p_envfit,
  width = 7,
  height = 5
)

# ============================================================
# 12. PCoA 可视化
# ============================================================

eig_value <- pcoa_arg$eig
eig_percent <- eig_value / sum(eig_value[eig_value > 0]) * 100

p_pcoa <- ggplot(pcoa_site, aes(PCoA1, PCoA2, color = ktype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = ktype), linetype = 2) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", round(eig_percent[1], 2), "%)"),
    y = paste0("PCoA2 (", round(eig_percent[2], 2), "%)"),
    color = "ARG group",
    title = "PCoA of overall ARG composition"
  )

ggsave(
  file.path(outp, "08_PCoA_ARG_profile_ktype.pdf"),
  p_pcoa,
  width = 6,
  height = 5
)

# ============================================================
# 13. 准备 ML 数据
# ============================================================

climate_dummy_ml <- model.matrix(~ climate_type, data = model_df)

if (ncol(climate_dummy_ml) > 1) {
  climate_dummy_ml <- climate_dummy_ml[, -1, drop = FALSE]
} else {
  climate_dummy_ml <- matrix(nrow = nrow(model_df), ncol = 0)
}

rownames(climate_dummy_ml) <- rownames(model_df)

ml_X <- cbind(
  model_df[, setdiff(colnames(model_df), "climate_type"), drop = FALSE],
  climate_dummy_ml
)

ml_X <- as.data.frame(ml_X)

ml_factor_map_num <- factor_group_map %>%
  filter(factor != "climate_type") %>%
  select(factor, original_factor, group)

ml_factor_map_climate <- tibble(
  factor = colnames(climate_dummy_ml),
  original_factor = "climate type",
  group = "Geo-climatic factors"
)

ml_factor_map <- bind_rows(
  ml_factor_map_num,
  ml_factor_map_climate
) %>%
  mutate(
    mlr_factor = paste0("X", row_number())
  )

colnames(ml_X) <- ml_factor_map$mlr_factor

ml_data <- data.frame(
  sample = rownames(arg_mat),
  PCoA1 = pcoa_arg$points[, 1],
  PCoA2 = pcoa_arg$points[, 2],
  ml_X,
  check.names = FALSE
)

write.csv(
  ml_factor_map,
  file.path(outp, "09_mlr_factor_group_map.csv"),
  row.names = FALSE
)

write.csv(
  ml_data,
  file.path(outp, "09_mlr_input_ARG_profile_PCoA.csv"),
  row.names = FALSE
)

# ============================================================
# 14. RF + SVM permutation importance
# ============================================================

rf_lrn <- makeLearner(
  "regr.randomForest",
  ntree = 1000,
  importance = TRUE
)

svm_lrn <- makeLearner(
  "regr.svm",
  kernel = "radial",
  cost = 1,
  gamma = 1 / length(ml_factor_map$mlr_factor),
  scale = TRUE
)

ml_importance_raw <- tibble()

n_perm <- 100

for (target_now in c("PCoA1", "PCoA2")) {
  
  task_data <- ml_data %>%
    select(all_of(target_now), all_of(ml_factor_map$mlr_factor))
  
  task_now <- makeRegrTask(
    id = paste0("ARG_profile_", target_now),
    data = task_data,
    target = target_now
  )
  
  set.seed(123)
  rf_model <- train(rf_lrn, task_now)
  
  set.seed(123)
  svm_model <- train(svm_lrn, task_now)
  
  for (model_now in c("RF", "SVM")) {
    
    if (model_now == "RF") {
      trained_model <- rf_model
    } else {
      trained_model <- svm_model
    }
    
    base_pred <- predict(trained_model, task = task_now)
    base_rmse <- performance(base_pred, measures = rmse)
    
    for (j in ml_factor_map$mlr_factor) {
      
      perm_values <- c()
      
      for (k in 1:n_perm) {
        
        tmp_data <- task_data
        tmp_data[[j]] <- sample(tmp_data[[j]])
        
        tmp_task <- makeRegrTask(
          id = paste0("tmp_", model_now, "_", target_now),
          data = tmp_data,
          target = target_now
        )
        
        tmp_pred <- predict(trained_model, task = tmp_task)
        tmp_rmse <- performance(tmp_pred, measures = rmse)
        
        perm_values <- c(perm_values, tmp_rmse - base_rmse)
      }
      
      ml_importance_raw <- bind_rows(
        ml_importance_raw,
        data.frame(
          target = target_now,
          model = model_now,
          mlr_factor = j,
          importance_mean = mean(perm_values),
          importance_sd = sd(perm_values),
          positive_rate = mean(perm_values > 0)
        )
      )
    }
  }
}

write.csv(
  ml_importance_raw,
  file.path(outp, "10_RF_SVM_raw_permutation_importance_ARG_profile.csv"),
  row.names = FALSE
)

# ============================================================
# 15. 按具体因子整合重要性
# ============================================================

ml_importance_factor <- ml_importance_raw %>%
  left_join(ml_factor_map, by = "mlr_factor") %>%
  group_by(target, model, original_factor, group) %>%
  summarise(
    importance_mean = sum(importance_mean, na.rm = TRUE),
    .groups = "drop"
  )

ml_importance_integrated <- ml_importance_factor %>%
  group_by(original_factor, group) %>%
  summarise(
    mean_importance = mean(importance_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rank_importance = rank(-mean_importance, ties.method = "average")
  ) %>%
  arrange(rank_importance)

write.csv(
  ml_importance_integrated,
  file.path(outp, "11_RF_SVM_integrated_importance_by_factor.csv"),
  row.names = FALSE
)

p_ml_factor <- ml_importance_integrated %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = mean_importance, fill = group)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Mean increase in RMSE after permutation",
    fill = "Factor group",
    title = "RF and SVM importance for overall ARG composition"
  )

ggsave(
  file.path(outp, "11_RF_SVM_integrated_importance_by_factor.pdf"),
  p_ml_factor,
  width = 7,
  height = 5
)

# ============================================================
# 16. 按因子组整合重要性
# ============================================================

ml_importance_group <- ml_importance_integrated %>%
  group_by(group) %>%
  summarise(
    total_importance = sum(mean_importance, na.rm = TRUE),
    mean_importance = mean(mean_importance, na.rm = TRUE),
    n_factor = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_importance))

write.csv(
  ml_importance_group,
  file.path(outp, "12_RF_SVM_integrated_importance_by_group.csv"),
  row.names = FALSE
)

p_ml_group <- ml_importance_group %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = group, y = total_importance, fill = group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "Total importance",
    title = "Grouped RF and SVM importance for overall ARG composition"
  )

ggsave(
  file.path(outp, "12_RF_SVM_integrated_importance_by_group.pdf"),
  p_ml_group,
  width = 6,
  height = 4
)




剔除一半为0的数据
# ============================================================
# 3. 整理整体 ARG subtype 丰度矩阵，并剔除一半以上样本为 0 的 subtype
# ============================================================

sample_cols <- intersect(colnames(nor_cell_sub_anno), rownames(factor_df))

arg_mat_raw <- nor_cell_sub_anno %>%
  select(subtype, all_of(sample_cols)) %>%
  group_by(subtype) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat_raw[is.na(arg_mat_raw)] <- 0

common_samples <- intersect(rownames(arg_mat_raw), rownames(factor_df))

arg_mat_raw <- arg_mat_raw[common_samples, , drop = FALSE]
factor_df <- factor_df[common_samples, , drop = FALSE]

# 统计每个 ARG subtype 在样本中的 0 值比例
arg_zero_stat <- tibble(
  subtype = colnames(arg_mat_raw),
  n_sample = nrow(arg_mat_raw),
  n_zero = colSums(arg_mat_raw == 0, na.rm = TRUE),
  zero_ratio = n_zero / n_sample,
  total_abundance = colSums(arg_mat_raw, na.rm = TRUE)
) %>%
  mutate(
    keep = zero_ratio <= 0.5 & total_abundance > 0
  ) %>%
  arrange(desc(zero_ratio))

write.csv(
  arg_zero_stat,
  file.path(outp, "01_ARG_subtype_zero_filter_summary.csv"),
  row.names = FALSE
)

# 保留在至少一半样本中出现的 ARG subtype
arg_mat <- arg_mat_raw[, arg_zero_stat$subtype[arg_zero_stat$keep], drop = FALSE]

write.csv(
  arg_mat %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_matrix_after_zero_filter.csv"),
  row.names = FALSE
)

# 输出筛选前后数量
filter_summary <- data.frame(
  raw_ARG_subtype_number = ncol(arg_mat_raw),
  retained_ARG_subtype_number = ncol(arg_mat),
  removed_ARG_subtype_number = ncol(arg_mat_raw) - ncol(arg_mat),
  zero_filter_threshold = "remove subtype with zero_ratio > 0.5"
)

write.csv(
  filter_summary,
  file.path(outp, "01_ARG_zero_filter_number_summary.csv"),
  row.names = FALSE
)

print(filter_summary)

# ============================================================
# 4. ARGs 标准化
#    arg_rel 用于 Bray-Curtis / PERMANOVA / PCoA / db-RDA
#    arg_hel 用于 varpart / RDA 类分析
# ============================================================

arg_rel <- decostand(arg_mat, method = "total")
arg_rel[is.na(arg_rel)] <- 0

arg_hel <- decostand(arg_rel, method = "hellinger")
arg_hel[is.na(arg_hel)] <- 0

arg_bray <- vegdist(arg_rel, method = "bray")

write.csv(
  arg_rel %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_relative_abundance_after_zero_filter.csv"),
  row.names = FALSE
)

write.csv(
  arg_hel %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_hellinger_after_zero_filter.csv"),
  row.names = FALSE
)



202600607补充 gdp
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(mlr)
library(randomForest)
library(e1071)

outp <- "outp/factor_ARGs_profile_grouped_GDP"
dir.create(outp, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 读取 factors0527.csv
# ============================================================

factors <- read.csv(
  "input/factors0527_lxc.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ============================================================
# 2. 因子分类
# ============================================================

environmental_vars <- c("As", "Hg", "P", "Cd", "Cr", "Pb", "N", "OM")

geo_climate_num_vars <- c(
  "Annual_average_temperature",
  "Annual_precipitation",
  "Annual_sunshine_hours",
  "longitude",
  "latitude"
)

geo_climate_cat_vars <- c("climate_type")

economic_vars <- c(
  "Green_area",
  "Per_capita_regional_GDP",
  "Total_population",
  "GDP"
)

factor_group_map_all <- tibble(
  factor = c(
    environmental_vars,
    geo_climate_num_vars,
    geo_climate_cat_vars,
    economic_vars
  ),
  original_factor = c(
    environmental_vars,
    "Annual average temperature",
    "Annual precipitation",
    "Annual sunshine hours",
    "longitude",
    "latitude",
    "climate type",
    "Green area",
    "Per capita regional GDP",
    "Total population",
    "GDP"
  ),
  group = c(
    rep("Environmental factors", length(environmental_vars)),
    rep("Geo-climatic factors", length(geo_climate_num_vars)),
    rep("Geo-climatic factors", length(geo_climate_cat_vars)),
    rep("Economic factors", length(economic_vars))
  )
)

# ============================================================
# 3. 整理 factors
# ============================================================

factor_df <- factors %>%
  mutate(across(everything(), ~na_if(as.character(.x), "/"))) %>%
  rename(
    climate_type = `climate type`,
    Annual_average_temperature = `Annual average temperature`,
    Annual_precipitation = `Annual precipitation`,
    Annual_sunshine_hours = `Annual sunshine hours`,
    Green_area = `Green area`,
    Per_capita_regional_GDP = `Per capita regional GDP`,
    Total_population = `Total population`
  ) %>%
  mutate(
    climate_type = trimws(climate_type)
  ) %>%
  mutate(across(
    all_of(c(environmental_vars, geo_climate_num_vars, economic_vars)),
    as.numeric
  )) %>%
  mutate(
    climate_type = replace_na(climate_type, "Unknown"),
    climate_type = factor(climate_type),
    ktype = factor(ktype)
  ) %>%
  column_to_rownames("sample")

# ============================================================
# 4. 整理 ARG subtype 矩阵，并剔除超过一半样本为 0 的 subtype
# ============================================================

sample_cols <- intersect(colnames(nor_cell_sub_anno), rownames(factor_df))

arg_mat_raw <- nor_cell_sub_anno %>%
  select(subtype, all_of(sample_cols)) %>%
  group_by(subtype) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat_raw[is.na(arg_mat_raw)] <- 0

common_samples <- intersect(rownames(arg_mat_raw), rownames(factor_df))

arg_mat_raw <- arg_mat_raw[common_samples, , drop = FALSE]
factor_df <- factor_df[common_samples, , drop = FALSE]

arg_zero_stat <- tibble(
  subtype = colnames(arg_mat_raw),
  n_sample = nrow(arg_mat_raw),
  n_zero = colSums(arg_mat_raw == 0, na.rm = TRUE),
  zero_ratio = n_zero / n_sample,
  total_abundance = colSums(arg_mat_raw, na.rm = TRUE)
) %>%
  mutate(
    keep = zero_ratio <= 0.9 & total_abundance > 0
  ) %>%
  arrange(desc(zero_ratio))

write.csv(
  arg_zero_stat,
  file.path(outp, "01_ARG_subtype_zero_filter_summary.csv"),
  row.names = FALSE
)

arg_mat <- arg_mat_raw[, arg_zero_stat$subtype[arg_zero_stat$keep], drop = FALSE]

filter_summary <- data.frame(
  raw_ARG_subtype_number = ncol(arg_mat_raw),
  retained_ARG_subtype_number = ncol(arg_mat),
  removed_ARG_subtype_number = ncol(arg_mat_raw) - ncol(arg_mat),
  zero_filter_threshold = "remove subtype with zero_ratio > 0.5"
)

write.csv(
  filter_summary,
  file.path(outp, "01_ARG_zero_filter_number_summary.csv"),
  row.names = FALSE
)

write.csv(
  arg_mat %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_matrix_after_zero_filter.csv"),
  row.names = FALSE
)

print(filter_summary)

# ============================================================
# 5. ARGs 标准化
# ============================================================

arg_rel <- decostand(arg_mat, method = "total")
arg_rel[is.na(arg_rel)] <- 0

arg_hel <- decostand(arg_rel, method = "hellinger")
arg_hel[is.na(arg_hel)] <- 0

arg_bray <- vegdist(arg_rel, method = "bray")

write.csv(
  arg_rel %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_relative_abundance_after_zero_filter.csv"),
  row.names = FALSE
)

write.csv(
  arg_hel %>% rownames_to_column("sample"),
  file.path(outp, "01_ARG_hellinger_after_zero_filter.csv"),
  row.names = FALSE
)
# ============================================================
# 6. 整理数值型因子
# ============================================================

all_num_vars <- c(environmental_vars, geo_climate_num_vars, economic_vars)

env_num_raw <- factor_df %>%
  select(any_of(all_num_vars))

num_keep <- names(which(colSums(!is.na(env_num_raw)) >= 0.7 * nrow(env_num_raw)))

env_num <- env_num_raw[, num_keep, drop = FALSE]

num_sd <- apply(env_num, 2, sd, na.rm = TRUE)
num_keep <- names(num_sd[num_sd > 0])

env_num <- env_num[, num_keep, drop = FALSE]

for (i in colnames(env_num)) {
  env_num[[i]][is.na(env_num[[i]])] <- median(env_num[[i]], na.rm = TRUE)
}

env_num_scale <- as.data.frame(scale(env_num))
rownames(env_num_scale) <- rownames(env_num)

geo_cat <- factor_df %>%
  select(any_of(geo_climate_cat_vars))

geo_cat$climate_type <- droplevels(factor(geo_cat$climate_type))

if (n_distinct(geo_cat$climate_type) > 1) {
  cat_keep <- "climate_type"
} else {
  cat_keep <- character(0)
}

model_df <- cbind(
  env_num_scale,
  geo_cat[, cat_keep, drop = FALSE]
)

rownames(model_df) <- rownames(factor_df)

factor_group_map <- factor_group_map_all %>%
  filter(factor %in% colnames(model_df))

write.csv(
  factor_group_map,
  file.path(outp, "01_factor_group_map_used.csv"),
  row.names = FALSE
)

write.csv(
  model_df %>% rownames_to_column("sample"),
  file.path(outp, "01_factor_matrix_scaled_used.csv"),
  row.names = FALSE
)
# ============================================================
# 7. 单因子 PERMANOVA
# ============================================================

permanova_single <- tibble()

for (i in colnames(model_df)) {
  
  tmp_env <- model_df[, i, drop = FALSE]
  colnames(tmp_env) <- "factor_value"
  
  if (is.numeric(tmp_env$factor_value)) {
    valid_factor <- sd(tmp_env$factor_value, na.rm = TRUE) > 0
  } else {
    valid_factor <- n_distinct(tmp_env$factor_value) > 1
  }
  
  if (valid_factor) {
    
    test <- adonis2(
      arg_bray ~ factor_value,
      data = tmp_env,
      permutations = 999
    )
    
    permanova_single <- bind_rows(
      permanova_single,
      data.frame(
        factor = i,
        R2 = test$R2[1],
        F_value = test$F[1],
        p_value = test$`Pr(>F)`[1]
      )
    )
  }
}

permanova_single <- permanova_single %>%
  left_join(factor_group_map, by = "factor") %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  permanova_single,
  file.path(outp, "02_single_factor_PERMANOVA_grouped.csv"),
  row.names = FALSE
)

p_single <- permanova_single %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = R2, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "PERMANOVA R2",
    fill = "Factor group",
    title = "Single-factor effects on overall ARG composition"
  )

ggsave(
  file.path(outp, "02_single_factor_PERMANOVA_grouped.pdf"),
  p_single,
  width = 7,
  height = 5
)
# ============================================================
# 8. 分组 PERMANOVA
# ============================================================

group_permanova <- tibble()

env_data <- model_df[, intersect(environmental_vars, colnames(model_df)), drop = FALSE]

if (ncol(env_data) > 0) {
  
  test_env <- adonis2(arg_bray ~ ., data = env_data, permutations = 999)
  
  group_permanova <- bind_rows(
    group_permanova,
    data.frame(
      group = "Environmental factors",
      n_factor = ncol(env_data),
      R2 = test_env$R2[1],
      F_value = test_env$F[1],
      p_value = test_env$`Pr(>F)`[1]
    )
  )
}

geo_data <- model_df[, intersect(c(geo_climate_num_vars, geo_climate_cat_vars), colnames(model_df)), drop = FALSE]

if (ncol(geo_data) > 0) {
  
  test_geo <- adonis2(arg_bray ~ ., data = geo_data, permutations = 999)
  
  group_permanova <- bind_rows(
    group_permanova,
    data.frame(
      group = "Geo-climatic factors",
      n_factor = ncol(geo_data),
      R2 = test_geo$R2[1],
      F_value = test_geo$F[1],
      p_value = test_geo$`Pr(>F)`[1]
    )
  )
}

econ_data <- model_df[, intersect(economic_vars, colnames(model_df)), drop = FALSE]

if (ncol(econ_data) > 0) {
  
  test_econ <- adonis2(arg_bray ~ ., data = econ_data, permutations = 999)
  
  group_permanova <- bind_rows(
    group_permanova,
    data.frame(
      group = "Economic factors",
      n_factor = ncol(econ_data),
      R2 = test_econ$R2[1],
      F_value = test_econ$F[1],
      p_value = test_econ$`Pr(>F)`[1]
    )
  )
}

group_permanova <- group_permanova %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(desc(R2))

write.csv(
  group_permanova,
  file.path(outp, "03_factor_group_PERMANOVA.csv"),
  row.names = FALSE
)

p_group <- group_permanova %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = group, y = R2, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = significance), hjust = -0.2, size = 5) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "PERMANOVA R2",
    title = "Grouped effects on overall ARG composition"
  )

ggsave(
  file.path(outp, "03_factor_group_PERMANOVA.pdf"),
  p_group,
  width = 6,
  height = 4
)
# ============================================================
# 9. db-RDA
# ============================================================

dbrda_full <- capscale(
  arg_bray ~ .,
  data = model_df
)

write.csv(
  as.data.frame(anova.cca(dbrda_full, permutations = 999)),
  file.path(outp, "04_dbRDA_overall_test.csv")
)

write.csv(
  as.data.frame(anova.cca(dbrda_full, by = "term", permutations = 999)),
  file.path(outp, "04_dbRDA_factor_test.csv")
)

site_score <- scores(dbrda_full, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(factor_df %>% rownames_to_column("sample"), by = "sample")

env_score <- scores(dbrda_full, display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column("factor") %>%
  left_join(factor_group_map, by = "factor")

dbrda_imp <- summary(dbrda_full)$cont$importance

dbrda_exp1 <- ifelse(ncol(dbrda_imp) >= 1, dbrda_imp[2, 1] * 100, 0)
dbrda_exp2 <- ifelse(ncol(dbrda_imp) >= 2, dbrda_imp[2, 2] * 100, 0)

p_dbrda <- ggplot(site_score, aes(CAP1, CAP2, color = ktype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = ktype), linetype = 2) +
  geom_segment(
    data = env_score,
    aes(
      x = 0,
      y = 0,
      xend = CAP1,
      yend = CAP2,
      linetype = group
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = env_score,
    aes(CAP1, CAP2, label = original_factor),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw() +
  labs(
    x = paste0("CAP1 (", round(dbrda_exp1, 2), "%)"),
    y = paste0("CAP2 (", round(dbrda_exp2, 2), "%)"),
    color = "ARG group",
    linetype = "Factor group",
    title = "db-RDA of overall ARG composition"
  )

ggsave(
  file.path(outp, "04_dbRDA_ARG_profile_grouped_factors.pdf"),
  p_dbrda,
  width = 7,
  height = 6
)

# ============================================================
# 10. PCoA + envfit
# ============================================================

pcoa_arg <- cmdscale(arg_bray, eig = TRUE, k = 2)

pcoa_site <- data.frame(
  sample = rownames(arg_mat),
  PCoA1 = pcoa_arg$points[, 1],
  PCoA2 = pcoa_arg$points[, 2]
) %>%
  left_join(factor_df %>% rownames_to_column("sample"), by = "sample")

envfit_res <- envfit(
  pcoa_arg$points,
  model_df,
  permutations = 999
)

envfit_vec <- as.data.frame(scores(envfit_res, display = "vectors")) %>%
  rownames_to_column("factor") %>%
  mutate(
    r2 = envfit_res$vectors$r,
    p_value = envfit_res$vectors$pvals,
    variable_type = "numeric"
  )

if (!is.null(envfit_res$factors)) {
  
  envfit_fac <- data.frame(
    factor = names(envfit_res$factors$r),
    r2 = as.numeric(envfit_res$factors$r),
    p_value = as.numeric(envfit_res$factors$pvals),
    variable_type = "categorical"
  )
  
} else {
  
  envfit_fac <- tibble(
    factor = character(),
    r2 = numeric(),
    p_value = numeric(),
    variable_type = character()
  )
}

envfit_df <- bind_rows(envfit_vec, envfit_fac) %>%
  left_join(factor_group_map, by = "factor") %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  envfit_df,
  file.path(outp, "05_envfit_ARG_profile_grouped_factors.csv"),
  row.names = FALSE
)

p_envfit <- envfit_df %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = r2, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "envfit r2",
    fill = "Factor group",
    title = "Environmental fitting for overall ARG composition"
  )

ggsave(
  file.path(outp, "05_envfit_ARG_profile_grouped_factors.pdf"),
  p_envfit,
  width = 7,
  height = 5
)

# ============================================================
# 11. PCoA 图
# ============================================================

eig_value <- pcoa_arg$eig
eig_percent <- eig_value / sum(eig_value[eig_value > 0]) * 100

p_pcoa <- ggplot(pcoa_site, aes(PCoA1, PCoA2, color = ktype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = ktype), linetype = 2) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", round(eig_percent[1], 2), "%)"),
    y = paste0("PCoA2 (", round(eig_percent[2], 2), "%)"),
    color = "ARG group",
    title = "PCoA of overall ARG composition"
  )

ggsave(
  file.path(outp, "06_PCoA_ARG_profile_ktype.pdf"),
  p_pcoa,
  width = 6,
  height = 5
)
# ============================================================
# 12. 计算筛选后 ARGs 总丰度
# ============================================================

arg_total <- data.frame(
  sample = rownames(arg_mat),
  total_ARG_abundance = rowSums(arg_mat, na.rm = TRUE)
)

min_pos <- min(arg_total$total_ARG_abundance[arg_total$total_ARG_abundance > 0], na.rm = TRUE)

arg_total <- arg_total %>%
  mutate(
    log10_total_ARG_abundance = log10(total_ARG_abundance + min_pos / 2)
  ) %>%
  left_join(
    factor_df %>% rownames_to_column("sample"),
    by = "sample"
  )

write.csv(
  arg_total,
  file.path(outp, "07_total_ARG_abundance_after_zero_filter.csv"),
  row.names = FALSE
)

# ============================================================
# 13. 多个因子与 ARGs 丰度的线性拟合
# ============================================================

linear_factor_df <- factor_df %>%
  select(any_of(all_num_vars))

linear_fit_summary <- tibble()

for (i in colnames(linear_factor_df)) {
  
  tmp <- data.frame(
    sample = rownames(linear_factor_df),
    ARG = arg_total$log10_total_ARG_abundance,
    factor_value = linear_factor_df[[i]]
  ) %>%
    drop_na()
  
  if (nrow(tmp) >= 5 && sd(tmp$ARG) > 0 && sd(tmp$factor_value) > 0) {
    
    fit <- lm(ARG ~ factor_value, data = tmp)
    fit_sum <- summary(fit)
    
    linear_fit_summary <- bind_rows(
      linear_fit_summary,
      data.frame(
        factor = i,
        n = nrow(tmp),
        slope = coef(fit)[2],
        intercept = coef(fit)[1],
        R2 = fit_sum$r.squared,
        adj_R2 = fit_sum$adj.r.squared,
        F_value = fit_sum$fstatistic[1],
        p_value = coef(fit_sum)[2, 4],
        AIC = AIC(fit)
      )
    )
  }
}

linear_fit_summary <- linear_fit_summary %>%
  left_join(factor_group_map_all, by = "factor") %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  linear_fit_summary,
  file.path(outp, "08_linear_fit_ARG_abundance_vs_factors.csv"),
  row.names = FALSE
)

# ============================================================
# 14. 线性拟合 R2 排序图
# ============================================================

p_linear_r2 <- linear_fit_summary %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = R2, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Linear model R2",
    fill = "Factor group",
    title = "Linear fitting between ARG abundance and factors"
  )

ggsave(
  file.path(outp, "08_linear_fit_ARG_abundance_R2.pdf"),
  p_linear_r2,
  width = 7,
  height = 5
)

# ============================================================
# 15. 多因子线性拟合散点图
# ============================================================

linear_plot_data <- arg_total %>%
  select(sample, ktype, log10_total_ARG_abundance) %>%
  left_join(
    factor_df %>%
      rownames_to_column("sample") %>%
      select(sample, any_of(all_num_vars)),
    by = "sample"
  ) %>%
  pivot_longer(
    cols = all_of(colnames(linear_factor_df)),
    names_to = "factor",
    values_to = "factor_value"
  ) %>%
  drop_na() %>%
  left_join(factor_group_map_all, by = "factor")

p_linear_all <- ggplot(
  linear_plot_data,
  aes(x = factor_value, y = log10_total_ARG_abundance)
) +
  geom_point(aes(color = ktype), size = 2.5, alpha = 0.85) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.7,
    color = "black"
  ) +
  facet_wrap(
    ~ original_factor,
    scales = "free_x",
    ncol = 4
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "log10(total ARG abundance)",
    color = "ARG group",
    title = "Linear fitting between ARG abundance and multiple factors"
  )

ggsave(
  file.path(outp, "09_linear_fit_ARG_abundance_multiple_factors.pdf"),
  p_linear_all,
  width = 12,
  height = 10
)
# ============================================================
# 16. 准备 ML 数据
# ============================================================

climate_dummy_ml <- model.matrix(~ climate_type, data = model_df)

if (ncol(climate_dummy_ml) > 1) {
  climate_dummy_ml <- climate_dummy_ml[, -1, drop = FALSE]
} else {
  climate_dummy_ml <- matrix(nrow = nrow(model_df), ncol = 0)
}

rownames(climate_dummy_ml) <- rownames(model_df)

ml_X <- cbind(
  model_df[, setdiff(colnames(model_df), "climate_type"), drop = FALSE],
  climate_dummy_ml
)

ml_X <- as.data.frame(ml_X)

ml_factor_map_num <- factor_group_map %>%
  filter(factor != "climate_type") %>%
  select(factor, original_factor, group)

ml_factor_map_climate <- tibble(
  factor = colnames(climate_dummy_ml),
  original_factor = "climate type",
  group = "Geo-climatic factors"
)

ml_factor_map <- bind_rows(
  ml_factor_map_num,
  ml_factor_map_climate
) %>%
  mutate(
    mlr_factor = paste0("X", row_number())
  )

colnames(ml_X) <- ml_factor_map$mlr_factor

ml_data <- data.frame(
  sample = rownames(arg_mat),
  PCoA1 = pcoa_arg$points[, 1],
  PCoA2 = pcoa_arg$points[, 2],
  ml_X,
  check.names = FALSE
)

write.csv(
  ml_factor_map,
  file.path(outp, "10_mlr_factor_group_map.csv"),
  row.names = FALSE
)

write.csv(
  ml_data,
  file.path(outp, "10_mlr_input_ARG_profile_PCoA.csv"),
  row.names = FALSE
)

# ============================================================
# 17. RF + SVM permutation importance
# ============================================================

rf_lrn <- makeLearner(
  "regr.randomForest",
  ntree = 1000,
  importance = TRUE
)

svm_lrn <- makeLearner(
  "regr.svm",
  kernel = "radial",
  cost = 1,
  gamma = 1 / length(ml_factor_map$mlr_factor),
  scale = TRUE
)

ml_importance_raw <- tibble()

n_perm <- 100

for (target_now in c("PCoA1", "PCoA2")) {
  
  task_data <- ml_data %>%
    select(all_of(target_now), all_of(ml_factor_map$mlr_factor))
  
  task_now <- makeRegrTask(
    id = paste0("ARG_profile_", target_now),
    data = task_data,
    target = target_now
  )
  
  set.seed(123)
  rf_model <- train(rf_lrn, task_now)
  
  set.seed(123)
  svm_model <- train(svm_lrn, task_now)
  
  for (model_now in c("RF", "SVM")) {
    
    if (model_now == "RF") {
      trained_model <- rf_model
    } else {
      trained_model <- svm_model
    }
    
    base_pred <- predict(trained_model, task = task_now)
    base_rmse <- performance(base_pred, measures = rmse)
    
    for (j in ml_factor_map$mlr_factor) {
      
      perm_values <- c()
      
      for (k in 1:n_perm) {
        
        tmp_data <- task_data
        tmp_data[[j]] <- sample(tmp_data[[j]])
        
        tmp_task <- makeRegrTask(
          id = paste0("tmp_", model_now, "_", target_now),
          data = tmp_data,
          target = target_now
        )
        
        tmp_pred <- predict(trained_model, task = tmp_task)
        tmp_rmse <- performance(tmp_pred, measures = rmse)
        
        perm_values <- c(perm_values, tmp_rmse - base_rmse)
      }
      
      ml_importance_raw <- bind_rows(
        ml_importance_raw,
        data.frame(
          target = target_now,
          model = model_now,
          mlr_factor = j,
          importance_mean = mean(perm_values),
          importance_sd = sd(perm_values),
          positive_rate = mean(perm_values > 0)
        )
      )
    }
  }
}

write.csv(
  ml_importance_raw,
  file.path(outp, "11_RF_SVM_raw_permutation_importance_ARG_profile.csv"),
  row.names = FALSE
)

ml_importance_factor <- ml_importance_raw %>%
  left_join(ml_factor_map, by = "mlr_factor") %>%
  group_by(target, model, original_factor, group) %>%
  summarise(
    importance_mean = sum(importance_mean, na.rm = TRUE),
    .groups = "drop"
  )

ml_importance_integrated <- ml_importance_factor %>%
  group_by(original_factor, group) %>%
  summarise(
    mean_importance = mean(importance_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rank_importance = rank(-mean_importance, ties.method = "average")
  ) %>%
  arrange(rank_importance)

write.csv(
  ml_importance_integrated,
  file.path(outp, "12_RF_SVM_integrated_importance_by_factor.csv"),
  row.names = FALSE
)

p_ml_factor <- ml_importance_integrated %>%
  mutate(original_factor = factor(original_factor, levels = rev(original_factor))) %>%
  ggplot(aes(x = original_factor, y = mean_importance, fill = group)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Mean increase in RMSE after permutation",
    fill = "Factor group",
    title = "RF and SVM importance for overall ARG composition"
  )

ggsave(
  file.path(outp, "12_RF_SVM_integrated_importance_by_factor.pdf"),
  p_ml_factor,
  width = 7,
  height = 5
)

ml_importance_group <- ml_importance_integrated %>%
  group_by(group) %>%
  summarise(
    total_importance = sum(mean_importance, na.rm = TRUE),
    mean_importance = mean(mean_importance, na.rm = TRUE),
    n_factor = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_importance))

write.csv(
  ml_importance_group,
  file.path(outp, "13_RF_SVM_integrated_importance_by_group.csv"),
  row.names = FALSE
)

p_ml_group <- ml_importance_group %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = group, y = total_importance, fill = group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "Total importance",
    title = "Grouped RF and SVM importance for overall ARG composition"
  )

ggsave(
  file.path(outp, "13_RF_SVM_integrated_importance_by_group.pdf"),
  p_ml_group,
  width = 6,
  height = 4
)

############################################################
## MetaCHIP + GTDB taxonomy + abundance
## HGT frequency framework
## Author: ChatGPT
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

############################################################
## 0. 参数区：只改这里
############################################################

hgt_file <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input/result/metachip/mc_20260414_112459_combined_pcofg_HGTs_ip90_al200bp_c75_ei80_f10kbp/mc_20260414_112459_pcofg_detected_HGTs.txt"
tax_file <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input/result/metachip/gtdb_taxonomy.tsv"

## 如果你已经有 MAG 丰度表，就填这里；没有就设为 NA
abun_file <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input/result/coverm/abundance.tsv"

outdir <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/result/metachip_hgt_summary"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## 主要分析层级：建议 genus
main_level <- "genus"

## 是否去掉“同层级内部”的 HGT（如 genus 内部）
## TRUE：只看不同 taxon 之间
## FALSE：保留同 taxon 内部事件，并单独用 n*(n-1) 或 choose(n,2) 做分母
exclude_within_taxon <- TRUE

## 高可信筛选条件
min_support_n <- 2
keep_end_match <- "no"
keep_full_length_match <- "no"

############################################################
## 1. 小工具函数
############################################################

parse_rank <- function(x, prefix) {
  pat <- paste0(prefix, "__[^;]*")
  out <- stringr::str_extract(x, pat)
  out <- sub(paste0("^", prefix, "__"), "", out)
  out[is.na(out) | out == ""] <- "Unknown"
  out
}

strip_support_pct <- function(x) {
  sub("\\(.*\\)$", "", x)
}

extract_genome_from_gene <- function(x) {
  sub("_[^_]+$", "", x)
}

calc_shannon_breadth <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  s <- sum(x)
  if (s <= 0) return(NA_real_)
  p <- x / s
  p <- p[p > 0]
  -sum(p * log(p))
}

calc_pair_tables <- function(hgt_df, tax_map, level_col = "genus", exclude_within = TRUE) {
  level_sym <- rlang::sym(level_col)
  
  ## 当前层级每个 taxon 含多少 genome
  taxon_size <- tax_map %>%
    distinct(genome, !!level_sym) %>%
    count(!!level_sym, name = "n_genomes") %>%
    rename(taxon = !!level_sym)
  
  dat <- hgt_df %>%
    transmute(
      donor_taxon = .data[[paste0("donor_", level_col)]],
      recipient_taxon = .data[[paste0("recipient_", level_col)]]
    ) %>%
    filter(!is.na(donor_taxon), !is.na(recipient_taxon))
  
  if (exclude_within) {
    dat <- dat %>% filter(donor_taxon != recipient_taxon)
  }
  
  ## directed
  directed_count <- dat %>%
    count(donor_taxon, recipient_taxon, name = "hgt_count") %>%
    left_join(taxon_size %>% rename(donor_taxon = taxon, n_donor = n_genomes), by = "donor_taxon") %>%
    left_join(taxon_size %>% rename(recipient_taxon = taxon, n_recipient = n_genomes), by = "recipient_taxon") %>%
    mutate(
      total_genome_pairs = if_else(
        donor_taxon == recipient_taxon,
        n_donor * (n_donor - 1),
        n_donor * n_recipient
      ),
      hgt_frequency = hgt_count / total_genome_pairs
    ) %>%
    arrange(desc(hgt_frequency), desc(hgt_count), donor_taxon, recipient_taxon)
  
  ## undirected
  undirected_count <- dat %>%
    mutate(
      taxon_a = if_else(donor_taxon <= recipient_taxon, donor_taxon, recipient_taxon),
      taxon_b = if_else(donor_taxon <= recipient_taxon, recipient_taxon, donor_taxon)
    ) %>%
    count(taxon_a, taxon_b, name = "hgt_count") %>%
    left_join(taxon_size %>% rename(taxon_a = taxon, n_a = n_genomes), by = "taxon_a") %>%
    left_join(taxon_size %>% rename(taxon_b = taxon, n_b = n_genomes), by = "taxon_b") %>%
    mutate(
      total_genome_pairs = if_else(
        taxon_a == taxon_b,
        choose(n_a, 2),
        n_a * n_b
      ),
      hgt_frequency = hgt_count / total_genome_pairs
    ) %>%
    arrange(desc(hgt_frequency), desc(hgt_count), taxon_a, taxon_b)
  
  list(
    directed = directed_count,
    undirected = undirected_count,
    taxon_size = taxon_size
  )
}

############################################################
## 2. 读取 MetaCHIP 结果
############################################################

hgt_raw <- readr::read_tsv(hgt_file, show_col_types = FALSE)

## 列名兼容
names(hgt_raw) <- make.names(names(hgt_raw))

## 解析 donor / recipient genome
hgt1 <- hgt_raw %>%
  mutate(
    gene1_genome = extract_genome_from_gene(Gene_1),
    gene2_genome = extract_genome_from_gene(Gene_2),
    direction_clean = strip_support_pct(direction),
    donor_genome = sub("-->.*$", "", direction_clean),
    recipient_genome = sub("^.*?-->", "", direction_clean),
    support_n = stringr::str_count(occurence.pcofg., "1"),
    support_p = substr(occurence.pcofg., 1, 1),
    support_c = substr(occurence.pcofg., 2, 2),
    support_o = substr(occurence.pcofg., 3, 3),
    support_f = substr(occurence.pcofg., 4, 4),
    support_g = substr(occurence.pcofg., 5, 5)
  )

## 高可信子集
hgt_use <- hgt1 %>%
  filter(
    support_n >= min_support_n,
    end_match == keep_end_match,
    full_length_match == keep_full_length_match
  )

readr::write_tsv(hgt1, file.path(outdir, "metachip_hgt_all.tsv"))
readr::write_tsv(hgt_use, file.path(outdir, "metachip_hgt_high_confidence.tsv"))

############################################################
## 3. 读取 GTDB taxonomy，并拆各层级
############################################################

tax_raw <- readr::read_tsv(
  tax_file,
  col_names = c("genome", "classification"),
  show_col_types = FALSE
)

tax_map <- tax_raw %>%
  mutate(
    domain  = parse_rank(classification, "d"),
    phylum  = parse_rank(classification, "p"),
    class   = parse_rank(classification, "c"),
    order   = parse_rank(classification, "o"),
    family  = parse_rank(classification, "f"),
    genus   = parse_rank(classification, "g"),
    species = parse_rank(classification, "s")
  )

readr::write_tsv(tax_map, file.path(outdir, "gtdb_taxonomy_parsed.tsv"))

############################################################
## 4. 给 HGT 表补 taxonomy
############################################################

donor_tax <- tax_map %>%
  rename_with(~ paste0("donor_", .x), -genome) %>%
  rename(donor_genome = genome)

recipient_tax <- tax_map %>%
  rename_with(~ paste0("recipient_", .x), -genome) %>%
  rename(recipient_genome = genome)

hgt_ann <- hgt_use %>%
  left_join(donor_tax, by = "donor_genome") %>%
  left_join(recipient_tax, by = "recipient_genome")

readr::write_tsv(hgt_ann, file.path(outdir, "metachip_hgt_high_confidence_annotated.tsv"))

############################################################
## 5. 各 taxonomic level 计算 HGT count / frequency
############################################################

levels_to_run <- c("phylum", "class", "order", "family", "genus")

pair_results <- purrr::map(levels_to_run, ~ calc_pair_tables(
  hgt_df = hgt_ann,
  tax_map = tax_map,
  level_col = .x,
  exclude_within = exclude_within_taxon
))
names(pair_results) <- levels_to_run

for (lv in levels_to_run) {
  readr::write_tsv(pair_results[[lv]]$taxon_size,
                   file.path(outdir, paste0(lv, "_taxon_size.tsv")))
  readr::write_tsv(pair_results[[lv]]$directed,
                   file.path(outdir, paste0(lv, "_pair_hgt_directed.tsv")))
  readr::write_tsv(pair_results[[lv]]$undirected,
                   file.path(outdir, paste0(lv, "_pair_hgt_undirected.tsv")))
}

############################################################
## 6. 主层级（默认 genus）做 donor / recipient hub 汇总
############################################################

main_directed <- pair_results[[main_level]]$directed

outgoing_summary <- main_directed %>%
  group_by(donor_taxon) %>%
  summarise(
    outgoing_hgt_count = sum(hgt_count, na.rm = TRUE),
    outgoing_hgt_frequency_sum = sum(hgt_frequency, na.rm = TRUE),
    n_recipient_taxa = n_distinct(recipient_taxon),
    .groups = "drop"
  ) %>%
  rename(taxon = donor_taxon)

incoming_summary <- main_directed %>%
  group_by(recipient_taxon) %>%
  summarise(
    incoming_hgt_count = sum(hgt_count, na.rm = TRUE),
    incoming_hgt_frequency_sum = sum(hgt_frequency, na.rm = TRUE),
    n_donor_taxa = n_distinct(donor_taxon),
    .groups = "drop"
  ) %>%
  rename(taxon = recipient_taxon)

taxon_hgt_summary <- full_join(outgoing_summary, incoming_summary, by = "taxon") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  arrange(desc(outgoing_hgt_count + incoming_hgt_count))

readr::write_tsv(taxon_hgt_summary,
                 file.path(outdir, paste0(main_level, "_hgt_summary.tsv")))

############################################################
## 7. 如果有 MAG 丰度表，计算 niche breadth
############################################################

if (!is.na(abun_file) && file.exists(abun_file)) {
  
  abun_raw <- readr::read_tsv(abun_file, show_col_types = FALSE)
  feature_col <- names(abun_raw)[1]
  
  ## MAG abundance table 转 long
  abun_long <- abun_raw %>%
    rename(feature = all_of(feature_col)) %>%
    filter(feature != "unmapped") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "abundance") %>%
    mutate(abundance = as.numeric(abundance))
  
  ## MAG level niche breadth
  mag_breadth <- abun_long %>%
    group_by(feature) %>%
    summarise(
      niche_breadth_shannon = calc_shannon_breadth(abundance),
      total_abundance = sum(abundance, na.rm = TRUE),
      n_samples_present = sum(abundance > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(genome = feature)
  
  readr::write_tsv(mag_breadth, file.path(outdir, "mag_niche_breadth.tsv"))
  
  ## genus level abundance
  genus_abun <- abun_long %>%
    rename(genome = feature) %>%
    left_join(tax_map %>% select(genome, genus), by = "genome") %>%
    filter(!is.na(genus), genus != "Unknown") %>%
    group_by(genus, sample) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
  
  genus_breadth <- genus_abun %>%
    group_by(genus) %>%
    summarise(
      niche_breadth_shannon = calc_shannon_breadth(abundance),
      total_abundance = sum(abundance, na.rm = TRUE),
      n_samples_present = sum(abundance > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(taxon = genus)
  
  readr::write_tsv(genus_breadth, file.path(outdir, "genus_niche_breadth.tsv"))
  
  ## genus niche breadth + HGT summary
  if (main_level == "genus") {
    genus_hgt_niche <- taxon_hgt_summary %>%
      left_join(genus_breadth, by = "taxon") %>%
      arrange(desc(outgoing_hgt_count + incoming_hgt_count))
    
    readr::write_tsv(genus_hgt_niche,
                     file.path(outdir, "genus_niche_breadth_vs_hgt.tsv"))
    
    ## 相关性
    cor_out <- suppressWarnings(
      cor.test(
        genus_hgt_niche$niche_breadth_shannon,
        genus_hgt_niche$outgoing_hgt_count,
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )
    
    cor_in <- suppressWarnings(
      cor.test(
        genus_hgt_niche$niche_breadth_shannon,
        genus_hgt_niche$incoming_hgt_count,
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )
    
    cor_tab <- tibble(
      comparison = c("niche_vs_outgoing_count", "niche_vs_incoming_count"),
      rho = c(unname(cor_out$estimate), unname(cor_in$estimate)),
      p_value = c(cor_out$p.value, cor_in$p.value)
    )
    readr::write_tsv(cor_tab, file.path(outdir, "genus_niche_hgt_correlation.tsv"))
    
    ## 图1：top outgoing genera
    p_out <- taxon_hgt_summary %>%
      slice_max(order_by = outgoing_hgt_count, n = 30) %>%
      mutate(taxon = forcats::fct_reorder(taxon, outgoing_hgt_count)) %>%
      ggplot(aes(x = taxon, y = outgoing_hgt_count)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(x = NULL, y = "Outgoing HGT count", title = "Top 30 donor genera")
    
    ggsave(file.path(outdir, "top30_outgoing_genera.pdf"),
           p_out, width = 8, height = 7)
    
    ## 图2：top incoming genera
    p_in <- taxon_hgt_summary %>%
      slice_max(order_by = incoming_hgt_count, n = 30) %>%
      mutate(taxon = forcats::fct_reorder(taxon, incoming_hgt_count)) %>%
      ggplot(aes(x = taxon, y = incoming_hgt_count)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(x = NULL, y = "Incoming HGT count", title = "Top 30 recipient genera")
    
    ggsave(file.path(outdir, "top30_incoming_genera.pdf"),
           p_in, width = 8, height = 7)
    
    ## 图3：niche breadth vs outgoing
    p_nb_out <- genus_hgt_niche %>%
      ggplot(aes(x = niche_breadth_shannon, y = outgoing_hgt_count)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw(base_size = 12) +
      labs(
        x = "Shannon niche breadth",
        y = "Outgoing HGT count",
        title = "Niche breadth vs outgoing HGT (genus)"
      )
    
    ggsave(file.path(outdir, "genus_niche_vs_outgoing_hgt.pdf"),
           p_nb_out, width = 6, height = 5)
    
    ## 图4：niche breadth vs incoming
    p_nb_in <- genus_hgt_niche %>%
      ggplot(aes(x = niche_breadth_shannon, y = incoming_hgt_count)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw(base_size = 12) +
      labs(
        x = "Shannon niche breadth",
        y = "Incoming HGT count",
        title = "Niche breadth vs incoming HGT (genus)"
      )
    
    ggsave(file.path(outdir, "genus_niche_vs_incoming_hgt.pdf"),
           p_nb_in, width = 6, height = 5)
  }
}

############################################################
## 8. 主层级 heatmap 输入表（top pairs）
############################################################

main_undir <- pair_results[[main_level]]$undirected

top_pairs <- main_undir %>%
  slice_max(order_by = hgt_frequency, n = 100)

readr::write_tsv(top_pairs,
                 file.path(outdir, paste0(main_level, "_top100_pair_hgt_frequency.tsv")))

############################################################
## 9. support_n 分布统计
############################################################

support_stat <- hgt1 %>%
  count(occurence.pcofg., support_n, end_match, full_length_match, name = "n_events") %>%
  arrange(desc(n_events))

readr::write_tsv(support_stat, file.path(outdir, "hgt_support_pattern_summary.tsv"))

############################################################
## 10. 运行信息
############################################################

message("All done.")
message("Output directory: ", outdir)
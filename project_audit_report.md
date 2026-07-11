# Project Audit Report

## 1. 项目结构概览

项目根目录下与分析直接相关的核心区域有：

- `input/`：原始输入、外部流程结果、人工整理表
- `output/`：当前较新的分析结果主目录
- `outp/`：较旧的分析结果主目录，仍被部分脚本继续消费
- 根目录若干 `*.R` / `*.r`：主分析脚本
- 根目录分析流程脚本：`readspipeline_cssd.sh`
- 根目录少量临时数据文件：如 `othersam2.csv`、`othersam3.csv`、`othersam4.csv`

共识别到 26 个 R 脚本，另有 1 个宏基因组主流程 shell 脚本。

## 2. R 脚本清单与大致输入/输出

### A. ARG 主线

- `数据整理.R`
  - 输入：`input/othersample*.csv`、`input/linshi*.csv/xlsx`、`input/sample.csv`
  - 输出：`input/othersam5.rda`
  - 作用：把分散样本表整理成后续很多脚本共用的基础对象

- `cssdarg.R`
  - 输入：`input/othersam5.rda`、`input/sarg/normalized_cell.*`、`input/sample.csv`
  - 输出：`outp/ARG_othersam5_3_剔除异常值/`
  - 作用：当前 ARG 丰度主结果之一；后续多个脚本直接读取该目录

- `ARG.R`
  - 输入：ARG abundance 与 sample 信息
  - 输出：`outp/arg_kmean`（脚本注释如此写），但现有结果更像沉淀在 `output/arg_kmean/`
  - 作用：ARG 聚类/可视化类分析

- `因子分析.R`
  - 输入：`input/factors0527_lxc.csv`、ARG abundance、sample metadata
  - 输出：
    - `outp/factor_ARGs_total_ML_0603`
    - `outp/factor_ARGs_profile`
    - `outp/factor_ARGs_profile_grouped`
    - `outp/factor_ARGs_profile_grouped_GDP`
  - 作用：环境/经济因子与 ARG 的统计和机器学习分析

### B. 细菌群落 / Kraken 主线

- `reads_kraken2.R`
  - 输入：`input/sample.csv`、`input/pluspf_taxid_7level_taxonomy.tsv`、`input/result/kraken2/bracken.all_levels.count.with_taxonomy.txt`
  - 输出：`output/read_kraken/`
  - 作用：构建 `microeco` 数据对象并计算多样性/网络基础结果

- `all_network.R`
  - 输入：`input/othersam5.rda`、Kraken 结果、`output/contig_taxid_tax_arg.rda`、病原菌/宿主分值相关结果
  - 输出：大量写入 `output/kraken_type1_distribution_network/`、`output/pathogen_*` 等
  - 作用：细菌网络、病原菌、宿主风险、部分整合分析的大脚本

- `all_net_arg_ana.R`
  - 输入：`input/othersam5.rda`、`outp/ARG_othersam5_3_剔除异常值/ARG_subtype_abundance_filtered_annotated.csv`、`output/kraken_type1_distribution_network/bacteria_species_network_type1/net.network.attribute.data.sample.csv`
  - 输出：`outp/ARG_network_mantel_type1/`
  - 作用：ARG 与网络属性 Mantel/相关分析

- `arg_net_mantel.R`
  - 输入：ARG abundance + network attributes
  - 输出：
    - `outp/arg_network_mantel/`
    - `outp/arg_network_mantel_by_ktype/`
    - `outp/arg_abundance_network/`
  - 作用：较早版本的 ARG-network Mantel 分析

- `Procrustes_analysis.r`
  - 输入：`input/othersam5.rda`、`outp/ARG_othersam5_3_剔除异常值/...`
  - 输出：`output/Procrustes_ARG_bacteria/`
  - 作用：ARG 与细菌群落 Procrustes 分析

### C. 病原菌主线

- `pathogen_env_analysis_type1.R`
  - 输入：
    - `input/othersam5.rda`
    - `input/pathogenic.csv`
    - `input/result/kraken2*/*.annotation.txt`
    - `input/result/kraken2*/*.count*.txt`
  - 输出：`output/result/pathogen_env_analysis_type1/`
  - 作用：不同环境中的病原菌丰度/丰富度分析
  - 特点：相对规范，支持项目根目录环境变量

- `water_sediment_ri_lefse.R`
  - 输入：病原菌/LEfSe/细菌对象等上游结果
  - 输出：
    - `output/result/lefse_rhizosphere_microbes/`
    - `output/result/rhizo_enriched_ARG_pathogen_summary/`
  - 作用：根际富集微生物与病原菌汇总

- `因子分析_病原菌_因子.R`
  - 输入：
    - `input/factors0527_lxc.csv`
    - `output/result/lefse_rhizosphere_microbes/`
    - `output/result/rhizo_enriched_ARG_pathogen_summary/`
    - `output/microeco_dataset_bacteria_type1.rds` 或同类对象
  - 输出：`output/result/pathogen_factor_fitting_factors0527/`
  - 作用：病原菌与环境因子拟合分析

- `rhizosphere_pathogen_sample_stacked_bar.R`
  - 输入：`input/factors0527_lxc.csv`、`input/othersam5.rda`、`input/pathogenic.csv`、Kraken 结果
  - 输出：`output/rhizosphere_pathogen_sample_stacked_bar/`
  - 作用：根际病原菌样本堆叠图

### D. contig / MAG / 宿主风险主线

- `contig_sarg.R`
  - 输入：
    - `input/result/kraken2/bracken.all_levels.count.with_taxonomy.txt`
    - `input/contig/ARGRANKER_DB.csv`
    - `input/mecha_col.rda`
    - `input/rank_col.rda`
  - 输出：
    - `output/contig_taxid_tax_arg.rda`
    - 以及 `output/pathogen_ARG_host_contig/` 相关结果
  - 作用：contig、病原菌、ARG、宿主关联
  - 特点：是多个下游脚本的重要中间层

- `MAG_ARG.R`
  - 输入：MAG / ARG / VF / MGE 相关结果
  - 输出：倾向于 `output/ARG_MGE_VF_host_score_*`
  - 作用：MAG-ARG 宿主风险整合

- `ARG_decrease_host_MGE_selection_pressure.R`
  - 输入：会联合读取 `input/` 与 `output/` 现有结果
  - 输出：`output/ARG_decrease_host_MGE_selection_pressure/`
  - 作用：把 ARG、宿主、MGE、选择压力整合成较新的结果链
  - 特点：结构较规范，包含 manifest 和 metadata 导出

- `rhizosphere_ARG_ktype_MAG_mechanism.R`
  - 输入：自身结果目录中的 CSV
  - 输出：`output/rhizosphere_ARG_ktype_MAG_mechanism/figures_effective_results/`
  - 作用：偏“结果可视化收尾”，不是原始计算入口

### E. 代谢组主线

- `metabolism_ana.R`
  - 输入：代谢物、宿主风险结果等
  - 输出：主要进入 `output/metabolite_microbe_host_score_analysis/`
  - 作用：代谢物与宿主/微生物关联分析的一部分

- `所有代谢与微生物的相关性分析.R`
  - 输入：
    - `input/sample.csv`
    - `input/metabolism.csv`
    - ARG abundance
    - `output/ARG_MGE_VF_host_score_strict/...`
    - `output/metabolite_microbe_host_score_analysis/...`
  - 输出：
    - `output/metabolite_pathogen_overall/`
    - `output/metabolite_pathogen_overall_optimized/`
    - `output/metabolite_pathogen_class_type/`
  - 作用：代谢物-病原菌/高风险宿主主线分析
  - 特点：明确存在“旧版结果 -> 优化版结果”的迭代关系

### F. vOTU / 病毒主线

- `vOTU.R`
  - 输入：`input/result/vOTU_formal/...`
  - 输出：
    - `output/vOTU_taxonomy_merged.tsv`
    - `output/vOTU_summary_annotated.tsv`
    - `output/vOTU_taxonomy_merge_stat.tsv`
    - 若干 vOTU lifestyle 顶层结果文件
  - 作用：整合 geNomad / CT3 / VIBRANT 的病毒分类结果

- `vibrant_amg.R`
  - 输入：
    - `input/result/salmon/gene.TPM`
    - `input/sample.csv`
    - `input/result/vOTU_formal/.../VIBRANT_AMG_individuals...tsv`
  - 输出：`output/vibrant_amg/`
  - 作用：AMG 注释、TPM、差异与富集分析

- `vOTU_host_arg_cytoscape.R`
  - 输入：
    - 自建目录 `vOTU_host_arg_cytoscape/node.csv`、`edge.csv`
    - `input/sarg/ARGRANKER_DB.csv`
    - `output/vOTU_taxonomy_merged.tsv`
    - `input/result/bin_MAG_function/MAG_gtdb_taxonomy.tsv`
  - 输出：写回 `vOTU_host_arg_cytoscape/`
  - 作用：Cytoscape 节点/边注释
  - 风险：输出位置不在 `output/` 下，而是项目根目录平级子目录

### G. 其他

- `map_arg.R`
  - 输入：依赖内存中的 `othersam5`、`nor_cell_sub_raw_all`，并使用 `input/NE2_HR_LC_SR_W_DR.tif`
  - 输出：`output/map_arg/`
  - 作用：ARG 地图绘制
  - 风险：不是完全独立脚本，似乎假定上游对象已存在

- `metachip.R`
  - 输入：`input/result/metachip/...`、`input/result/coverm/abundance.tsv`
  - 输出：`output/result/metachip_hgt_summary/`
  - 作用：HGT 结果汇总

- `hu_huanyong_line_rhizosphere_comparison.R`
  - 输入：`input/othersam5.rda`、`input/sample.csv`、若干结果对象
  - 输出：`output/hu_huanyong_line_rhizosphere_comparison/`
  - 作用：胡焕庸线比较分析

- `qinling_huaihe_line_rhizosphere_comparison.R`
  - 输入/输出：与 `hu_huanyong_line_rhizosphere_comparison.R` 高度耦合
  - 输出：`output/qinling_huaihe_line_rhizosphere_comparison/`
  - 作用：更像补充或派生比较脚本

## 3. 宏基因组主流程脚本 `readspipeline_cssd.sh`

### 3.1 脚本定位

- 文件：`readspipeline_cssd.sh`
- 性质：当前 CSSD 宏基因组分析的上游主流程脚本
- 作用：把原始双端 reads 从预处理一路推进到物种注释、丰度估计、组装、基因定量、CoverM、MetaCHIP、vOTU 等结果层
- 与本项目关系：`input/result/` 下的很多目录和文件都明显对应此脚本的产物，而不是手工生成

### 3.2 主要流程阶段

- 数据预处理
  - 原始数据从外部目录复制到工作目录
  - 批量重命名 `*.fq.gz` 为标准化的 `*_1.fastq.gz` / `*_2.fastq.gz`
  - 根据 `seq/` 文件名生成 `result/metadata.txt`
  - 用 `seqkit stat` 记录原始测序文件统计

- 质量控制
  - 使用 `fastp`
  - 输出 `temp/qc/*.json`、`*.html`、`*.log`
  - 生成 `result/qc/fastp.txt`
  - 将质控后的 `*.fastq` 转移到 `temp/hr/`

- Kraken2 物种注释
  - 输入：`temp/hr/{sample}_?.fastq`
  - 参考库：`${db}/kraken2/${type}`
  - 输出：
    - `temp/kraken2/{sample}.report`
    - `temp/kraken2/{sample}.output`
    - `temp/kraken2/{sample}.mpa`
    - `result/kraken2/tax_count.mpa`

- Bracken 丰度估计
  - 对 `D/P/C/O/F/G/S` 分类水平循环计算
  - 输出：
    - `temp/bracken/{sample}_{tax}.brk`
    - `temp/bracken/{sample}_{tax}.report`
    - `result/kraken2/bracken.*.count.txt`
    - `result/kraken2/bracken.*.relative.txt`
    - `result/kraken2/bracken.all_levels.*`

- Bracken 过滤与注释补回
  - 调用 `filter_feature_table_tang.R`
  - 调用 `filter_feature_table.R`
  - 调用 `add_bracken_taxonomy_tang.R`
  - 输出：
    - `bracken.${tax}.count.${prop}.txt`
    - `bracken.${tax}.relative.${prop}.txt`
    - `bracken.${tax}.count.${prop}.with_taxonomy.txt`
    - `bracken.${tax}.relative.${prop}.with_taxonomy.txt`

- 多样性与可视化
  - `alpha_diversity.py`
  - `beta_diversity.py`
  - `kreport2krona.py` + `ktImportText`
  - 输出：
    - `result/kraken2/alpha*.txt`
    - `result/kraken2/beta*.txt`
    - `result/kraken2/krona.*.html`

- 组装 / binning / 后续分析
  - 脚本后半段继续推进到组装、bin、vOTU、`salmon`、`coverm`、`MetaCHIP`
  - 明显对应项目中的：
    - `input/result/salmon/`
    - `input/result/coverm/`
    - `input/result/metachip/`
    - `input/result/bin_*`
    - `input/result/vOTU_formal/`

### 3.3 与当前项目目录的映射关系

从脚本内容看，下面这些项目输入目录大概率来自该 pipeline 或其直接派生流程：

- `input/result/kraken2/`
- `input/result/kraken2_106/`
- `input/result/kraken2_198/`
- `input/result/kraken2_ld300/`
- `input/result/kraken2_ldnc/`
- `input/result/salmon/`
- `input/result/coverm/`
- `input/result/metachip/`
- `input/result/bin_MAG_function/`
- `input/result/bin_MGE/`
- `input/result/bin_SARG/`
- `input/result/bin_VFDB/`
- `input/result/bin_intersect/`
- `input/result/vOTU_formal/`

这意味着本仓库中的很多 `input/result/*` 并不是“原始输入”，而是上游流程沉淀下来的中间或终端结果。

### 3.4 发现的问题

- 脚本中文注释本身存在明显乱码，说明文件编码已受损
- 脚本不是单一清晰版本，而是长期追加形成，存在多个历史块并存
- Bracken 段至少有两套实现并存
  - 一套是较新的 `20260430` 改进版
  - 后面还有较旧的 Bracken 汇总/过滤/多样性段
- 变量反复重定义
  - 如 `wd`、`cpu`、`db`、`output_dir`、`prop`
- 有不少硬编码服务器路径
  - 例如 `/home/tang/cssd`
  - 例如 `/home/tang/db`
  - 例如 `/EasyMicrobiome/...`
- 混合了“正式流程”和“调试/补跑记录”
  - 例如 `metadata1.txt`、`metadata2.txt`、`metadata4.txt`
  - 例如 `20260227 cpu02 ...`
- 同一脚本同时覆盖
  - 读长分类
  - 组装
  - vOTU
  - `salmon`
  - `coverm`
  - `MetaCHIP`
  - 可维护性较差，复跑风险较高

### 3.5 对审计的意义

- 该脚本应视为“项目的上游生产线说明书”
- 它解释了为什么 `input/result/` 下会同时存在 `kraken2`、`salmon`、`coverm`、`metachip`、`vOTU_formal` 等多条结果链
- 后续若要清理 `input/result/` 中的旧结果、重复结果或未调用文件，不能只看 R 脚本，还应结合这个 pipeline 判断哪些是上游必需产物

## 4. input 目录输入文件概况

### 关键基础表

- `input/sample.csv`
- `input/othersam5.rda`
- `input/factors0527_lxc.csv`
- `input/pathogenic.csv`
- `input/metabolism.csv`
- `input/pluspf_taxid_7level_taxonomy.tsv`

### 关键子目录

- `input/result/`
  - Kraken/Bracken
  - MAG/function/MGE/VFDB/SARG
  - metachip
  - vOTU_formal
  - coverm / salmon
- `input/contig/`
- `input/sarg/`
- `input/vibrant/`

### 疑似人工整理或临时输入

- `input/linshi.csv`
- `input/linshi.xlsx`
- `input/linshi2.csv`
- `input/linshi2.xlsx`
- `input/新建 Microsoft Excel 工作表.xlsx`
- `input/病原菌分类整理.xlsx`

这些更像手工过渡文件，不像稳定主输入。

## 5. output / outp 结果目录关系

### 当前较新的主结果目录

主要集中在 `output/`：

- `kraken_type1_distribution_network`
- `pathogen_distribution_type1`
- `pathogen_ARG_host_contig`
- `ARG_decrease_host_MGE_selection_pressure`
- `metabolite_microbe_host_score_analysis`
- `metabolite_pathogen_overall_optimized`
- `metabolite_pathogen_class_type`
- `vibrant_amg`
- `map_arg`
- `read_kraken`
- `result/pathogen_env_analysis_type1`
- `result/lefse_rhizosphere_microbes`
- `result/pathogen_factor_fitting_factors0527`

### 较旧或历史结果目录

主要集中在 `outp/`：

- `ARG_othersam5`
- `ARG_othersam5_1`
- `ARG_othersam5_2`
- `ARG_othersam5_3_剔除异常值`
- `arg_network_mantel_by_ktype`
- `arg_abundance_network`
- `factor_ARGs*`

其中当前脚本实际仍在读取的旧目录，最重要的是：

- `outp/ARG_othersam5_3_剔除异常值/`
- `outp/ARG_network_mantel_type1/`
- `outp/factor_ARGs*`

## 6. 旧版 / 新版 / 重复 / 疑似未调用结果

### 明确存在版本迭代关系

- `output/metabolite_pathogen_overall/`
  - 被 `所有代谢与微生物的相关性分析.R` 后半段继续读取
  - 明显是旧版整体分析

- `output/metabolite_pathogen_overall_optimized/`
  - 明显是新版优化结果
  - 当前更值得保留

- `outp/ARG_othersam5/`、`outp/ARG_othersam5_1/`、`outp/ARG_othersam5_2/`、`outp/ARG_othersam5_3_剔除异常值/`
  - 当前被后续脚本消费的是 `_3_剔除异常值`
  - `_1`、`_2` 很像中间迭代版本

### 疑似重复或历史残留

- `output/ARG_host_score_gg/`
- `output/ARG_MGE_VF_host_score_gg/`
- `output/sarg备份/`

在脚本全文搜索中没有发现这些目录被当前 R 脚本引用，属于“疑似孤立结果”。

### 顶层散落结果

`output/` 根目录直接存在一批文件，而不是落在子目录里：

- `vOTU_*`
- `Top100_species_*`
- `contig_taxid_tax_arg.rda`

其中：

- `contig_taxid_tax_arg.rda` 是重要中间文件，确实被脚本读取
- `vOTU_*` 是 `vOTU.R` 主输出，合理
- `Top100_species_*` 更像早期绘图产物，当前未见明确下游引用

### 明显异常文件

- `output/Top100_species_heatmap_ComplexHeatmap_Host_ktype.png` 大小为 0
  - 属于失败或未完整生成的结果，建议标记

## 7. 额外风险点

- 多个脚本曾长期依赖绝对路径，虽然目前已统一到当前项目路径，但仍不如相对路径稳妥
- 新旧结果目录混用
  - 新脚本会读取旧 `outp/` 结果，说明暂时不能简单删除 `outp/`
- `readspipeline_cssd.sh` 自身存在乱码、历史块叠加和重复实现
- `vOTU_host_arg_cytoscape.R` 的输出路径仍不标准
- `map_arg.R` 对内存对象的隐式依赖较强

## 8. 当前建议结论

- 现行主链建议优先围绕 `input/ -> output/` 梳理
- `input/result/` 应视为 pipeline 产生的上游结果层，而不是单纯原始输入
- `outp/` 应视为“历史但仍部分被依赖”的旧结果层，不能直接视为废弃
- 最值得优先人工复核的疑似旧版/重复目录：
  - `outp/ARG_othersam5_1`
  - `outp/ARG_othersam5_2`
  - `output/ARG_host_score_gg`
  - `output/ARG_MGE_VF_host_score_gg`
  - `output/sarg备份`
- 最值得优先修复的脚本问题：
  - `readspipeline_cssd.sh` 的编码与重复历史块
  - 将绝对路径进一步改成项目相对路径或 `project_root`
  - `vOTU_host_arg_cytoscape.R` 的非标准输出位置
  - `map_arg.R` 对内存对象的隐式依赖


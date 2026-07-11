
# 一、数据预处理 Data preprocessing

## 1.1 准备工作 Preparing
数据位置
/home/data/project/bn251213-lxh3/BN251119BJ01S27N5-83/30G-30samples
# 复制文件夹（-r 表示递归复制，用于文件夹及内部所有内容）
cp -r /home/data/project/bn251213-lxh3/BN251119BJ01S27N5-83/30G-30samples/rawdata ./cssd/
#将所有文件移动到seq
find ./rawdata -type f -name "*.fq.gz" -exec mv {} /home/tang/cssd/seq/ \;


    wd=/home/tang/cssd
	cpu=76
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
    mkdir -p seq temp result
 	# 添加分析所需的软件、脚本至环境变量，添加至~/.bashrc中自动加载，其中EasyMicrobiome文件夹放在根目录！！！
    PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/EasyMicrobiome/linux:/EasyMicrobiome/script
    echo $PATH
#修改文件名 
	cd ${wd}/seq
	# 切换到目标目录（确保当前在 seq 文件夹下）
cd /home/tang/cssd/seq

# 批量修改文件名（分两步替换，避免冲突）
for file in *.fq.gz; do
  # 第一步：替换 R1/R2 为 1/2；第二步：替换 fq.gz 为 fastq.gz
  new_name=$(echo "$file" | sed 's/R1/1/g; s/R2/2/g; s/fq.gz/fastq.gz/g')
  # 执行重命名（-n 参数用于测试，确认无误后删除 -n 再执行）
  mv -n "$file" "$new_name"
done

	#clean.sample.R1.fastq.gz 会被重命名为 sample_1.fastq.gz
	#clean.sample.R2.fastq.gz 会被重命名为 sample_2.fastq.gz



    #(c) 上传元数据metadata.txt至result目录，或者基于seq中文件名称快速生成metadata.txt
	cd ${wd}
	# 快速读取文件生成样本ID列表再继续编写
    ls $wd/seq/ | grep _1 | cut -f 1 -d '_' | sort -u | sed '1 i SampleID' > result/metadata.txt
    # 预览
    cat result/metadata.txt


    # 元数据细节优化
    # 转换Windows回车为Linux换行，去除空格
    sed -i 's/\r//;s/ //g' result/metadata.txt
    cat -A result/metadata.txt

    #(d) ls查看文件大小，-l 列出详细信息 (l: list)，-sh 显示人类可读方式文件大小 (s: size; h: human readable)
    ls -lsh seq/*.fastq.gz
    # 统计
    time seqkit stat seq/*.fastq.gz > result/seqkit.txt


	#**工作目录和文件结构总结**
    # ├── result
    # │   └── metadata.txt
    # ├── seq
    # │   ├── C1_1.fastq.gz
    # │   ├── ...
    # │   └── C1_2.fastq.gz
    # └── temp

## 1.2 Fastp质量控制 Quality Control

    # 创建目录，记录软件版本和引文
	cd $wd
    mkdir -p temp/qc result/qc
    # 多样本并行，此步占用原始数据5x空间
    # -j 2: 表示同时处理2个样本；j3,18s,8m
time tail -n+2 result/metadata.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j 30 \
  "fastp -i seq/{1}_1.fastq.gz -I seq/{1}_2.fastq.gz \
    -j temp/qc/{1}_fastp.json -h temp/qc/{1}_fastp.html \
    -o temp/qc/{1}_1.fastq  -O temp/qc/{1}_2.fastq \
    > temp/qc/{1}.log 2>&1"

    # 质控后结果汇总
    echo -e "SampleID\tRaw\tClean" > temp/fastp
    for i in `tail -n+2 result/metadata.txt|cut -f1`;do
        echo -e -n "$i\t" >> temp/fastp
        grep 'total reads' temp/qc/${i}.log|uniq|cut -f2 -d ':'|tr '\n' '\t' >> temp/fastp
        echo "" >> temp/fastp
        done
    sed -i 's/ //g;s/\t$//' temp/fastp
    # 按metadata排序
    awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' temp/fastp result/metadata.txt \
      > result/qc/fastp.txt
    cat result/qc/fastp.txt
    
	#跳过后直接将temp/qc/*fastq文件剪切至temp/hr
	mkdir -p temp/hr
	mv  $wd/temp/qc/*.fastq $wd/temp/hr
# 二、基于读长分析 Read-based (HUMAnN3+MetaPhlAn4+Kraken2)
    


## 2.7 Kraken2+Bracken物种注释和丰度估计
#1224 cpu01
#Kraken2可以快速完成读长(read)层面的物种注释和定量，还可以进行重叠群(contig)、基因(gene)、宏基因组组装基因组(MAG/bin)层面的序列物种注释。
tmux new -s cssdkraken
 wd=/home/tang/cssd
	cpu=76
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
    mkdir -p seq temp result
 	# 添加分析所需的软件、脚本至环境变量，添加至~/.bashrc中自动加载，其中EasyMicrobiome文件夹放在根目录！！！
    PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/EasyMicrobiome/linux:/EasyMicrobiome/script
    echo $PATH
		cd $wd
    # 启动kraken2工作环境
    conda activate kraken2
    # 记录软件版本
    kraken2 --version # 2.1.2
    mkdir -p temp/kraken2

### Kraken2物种注释

# 输入：temp/qc/{1}_?.fastq 质控后的数据，{1}代表样本名；
# 参考数据库：-db ${db}/kraken2/pluspf16g/
# 输出结果：每个样本单独输出，temp/kraken2/中的{1}_report和{1}_output
# 整合物种丰度表输出结果：result/kraken2/taxonomy_count.txt 

#多样本并行生成report，1样本8线程逐个运行，内存大但速度快，不建议用多任务并行
	type=pluspf
    for i in `tail -n+2 result/metadata.txt | cut -f1`;do
      kraken2 --db ${db}/kraken2/${type} \
      --paired temp/hr/${i}_?.fastq \
      --threads ${cpu} --use-names --report-zero-counts \
      --report temp/kraken2/${i}.report \
      --output temp/kraken2/${i}.output; done

     
#使用krakentools转换report为mpa格式

    for i in `tail -n+2 result/metadata.txt | cut -f1`;do
      kreport2mpa.py -r temp/kraken2/${i}.report \
        --display-header -o temp/kraken2/${i}.mpa; done

#合并样本为表格

    mkdir -p result/kraken2
    # 输出结果行数相同，但不一定顺序一致，要重新排序
    tail -n+2 result/metadata.txt | cut -f1 | rush -j 1 \
      'tail -n+2 temp/kraken2/{1}.mpa | LC_ALL=C sort | cut -f 2 | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count '
    # 提取第一样本品行名为表行名
    header=`tail -n 1 result/metadata.txt | cut -f 1`
    echo $header
    tail -n+2 temp/kraken2/${header}.mpa | LC_ALL=C sort | cut -f 1 | \
      sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
    head -n3 temp/kraken2/0header_count
    # paste合并样本为表格
    ls temp/kraken2/*count
    paste temp/kraken2/*count > result/kraken2/tax_count.mpa
    # 检查表格及统计
    csvtk -t stat result/kraken2/tax_count.mpa
    head -n 5 result/kraken2/tax_count.mpa
#1224 cpu01 end——————————————————————————————————

#bracken改进20260430————————————

# ============================================================
# 4. Bracken丰度估计
# ============================================================
# Bracken分类水平：
#   D = Domain
#   P = Phylum
#   C = Class
#   O = Order
#   F = Family
#   G = Genus
#   S = Species
#
# 重要：
#   bracken -t 不是线程数，而是最小reads阈值。
#   不要再写 -t ${cpu}。
# ============================================================

mkdir -p temp/bracken result/kraken2

readLen=150
prop=0.2
type=pluspf
threshold=0

for tax in D P C O F G S; do

  for i in $(tail -n+2 result/metadata.txt | cut -f1); do

    echo "Running Bracken: sample=${i}, level=${tax}"

    bracken \
      -d "${db}/kraken2/${type}/" \
      -i "temp/kraken2/${i}.report" \
      -r "${readLen}" \
      -l "${tax}" \
      -t "${threshold}" \
      -o "temp/bracken/${i}_${tax}.brk" \
      -w "temp/bracken/${i}_${tax}.report"

    if [ $? -ne 0 ]; then
      echo "ERROR: Bracken failed for sample ${i} at tax level ${tax}" >&2
    fi

  done

done


# ============================================================
# 5. 汇总Bracken结果
# ============================================================

db=/home/tang/db
input_dir="temp/bracken"
output_dir="result/kraken2"
metadata="result/metadata.txt"

mkdir -p "${output_dir}"

Rscript "$db/EasyMicrobiome/script/filter_feature_table_tang.R" \
  "${input_dir}" \
  "${output_dir}" \
  "${metadata}" \
  "D,P,C,O,F,G,S"

# ============================================================
# 6. 检查Bracken汇总结果
# ============================================================

echo "Check Bracken output tables"

csvtk -t stat result/kraken2/bracken.all_levels.long.txt
csvtk -t stat result/kraken2/bracken.all_levels.count.txt
csvtk -t stat result/kraken2/bracken.all_levels.relative.txt

head -n 5 result/kraken2/bracken.all_levels.long.txt
head -n 5 result/kraken2/bracken.all_levels.count.with_taxonomy.txt
head -n 5 result/kraken2/bracken.all_levels.relative.with_taxonomy.txt

# 查看每个分类水平的分类单元数量
cut -f3 result/kraken2/bracken.all_levels.long.txt | sort | uniq -c

# ============================================================
# 7. 过滤Bracken丰度表，并补回Taxonomy注释
# ============================================================
# 输入：
#   bracken.${tax}.count.txt
#   bracken.${tax}.relative.txt
#   bracken.all_levels.count.txt
#   bracken.all_levels.relative.txt
#
# 输出：
#   1）纯矩阵过滤结果：
#      bracken.${tax}.count.${prop}.txt
#      bracken.${tax}.relative.${prop}.txt
#
#   2）带Taxonomy注释的过滤结果：
#      bracken.${tax}.count.${prop}.with_taxonomy.txt
#      bracken.${tax}.relative.${prop}.with_taxonomy.txt
#
# 注意：
#   filter_feature_table.R 要求第一列为FeatureID，后面为样本列；
#   因此先过滤纯矩阵，再根据 annotation 文件补回分类注释。
# ============================================================

cd /home/tang/cssd

db=/home/tang/db
output_dir="result/kraken2"
prop=0.2

filter_script="${db}/EasyMicrobiome/script/filter_feature_table.R"
add_anno_script="${db}/EasyMicrobiome/script/add_bracken_taxonomy_tang.R"

mkdir -p "${output_dir}"


# ------------------------------------------------------------
# 7.1 定义过滤 + 补注释函数
# ------------------------------------------------------------

filter_and_add_taxonomy() {

  input_file="$1"
  annotation_file="$2"
  output_file="$3"
  output_tax_file="$4"
  label="$5"

  echo "--------------------------------------------------"
  echo "Filtering ${label}"
  echo "Input: ${input_file}"

  if [ ! -s "${input_file}" ]; then
    echo "WARNING: input file not found or empty: ${input_file}" >&2
    return 0
  fi

  if [ ! -s "${annotation_file}" ]; then
    echo "WARNING: annotation file not found or empty: ${annotation_file}" >&2
    return 0
  fi

  Rscript "${filter_script}" \
    -i "${input_file}" \
    -p "${prop}" \
    -o "${output_file}"

  if [ ! -s "${output_file}" ]; then
    echo "WARNING: filtered output not found or empty: ${output_file}" >&2
    return 0
  fi

  Rscript "${add_anno_script}" \
    "${output_file}" \
    "${annotation_file}" \
    "${output_tax_file}"

  echo "Finished: ${label}"
}


# ------------------------------------------------------------
# 7.2 过滤每个分类水平 D/P/C/O/F/G/S
# ------------------------------------------------------------

for tax in D P C O F G S; do

  annotation_file="${output_dir}/bracken.${tax}.annotation.txt"

  filter_and_add_taxonomy \
    "${output_dir}/bracken.${tax}.count.txt" \
    "${annotation_file}" \
    "${output_dir}/bracken.${tax}.count.${prop}.txt" \
    "${output_dir}/bracken.${tax}.count.${prop}.with_taxonomy.txt" \
    "Bracken ${tax} count table"

  filter_and_add_taxonomy \
    "${output_dir}/bracken.${tax}.relative.txt" \
    "${annotation_file}" \
    "${output_dir}/bracken.${tax}.relative.${prop}.txt" \
    "${output_dir}/bracken.${tax}.relative.${prop}.with_taxonomy.txt" \
    "Bracken ${tax} relative table"

done


# ------------------------------------------------------------
# 7.3 过滤所有分类水平合并表
# ------------------------------------------------------------
# 注意：
#   all_levels 表同时包含 D/P/C/O/F/G/S，
#   不能把其每个样本列加和理解为 1 或 100%。
#   该表适合查询、筛选、注释整合；
#   不建议直接用于 alpha/beta 多样性分析。
# ------------------------------------------------------------

annotation_file="${output_dir}/bracken.all_levels.annotation.txt"

filter_and_add_taxonomy \
  "${output_dir}/bracken.all_levels.count.txt" \
  "${annotation_file}" \
  "${output_dir}/bracken.all_levels.count.${prop}.txt" \
  "${output_dir}/bracken.all_levels.count.${prop}.with_taxonomy.txt" \
  "Bracken all levels count table"

filter_and_add_taxonomy \
  "${output_dir}/bracken.all_levels.relative.txt" \
  "${annotation_file}" \
  "${output_dir}/bracken.all_levels.relative.${prop}.txt" \
  "${output_dir}/bracken.all_levels.relative.${prop}.with_taxonomy.txt" \
  "Bracken all levels relative table"


# ------------------------------------------------------------
# 7.4 检查输出结果
# ------------------------------------------------------------

echo "=================================================="
echo "Filtered Bracken output files:"
echo "=================================================="

ls -lh "${output_dir}"/bracken.*.${prop}.txt
ls -lh "${output_dir}"/bracken.*.${prop}.with_taxonomy.txt

echo "=================================================="
echo "Check filtered all-levels table:"
echo "=================================================="

csvtk -t stat "${output_dir}/bracken.all_levels.relative.${prop}.txt"
csvtk -t stat "${output_dir}/bracken.all_levels.relative.${prop}.with_taxonomy.txt"

head -n 5 "${output_dir}/bracken.all_levels.relative.${prop}.with_taxonomy.txt"

echo "Bracken filtering and taxonomy annotation finished."

# ============================================================
# 8. 常用筛选示例
# ============================================================

# 从long总表中提取属水平
awk -F '\t' 'NR==1 || $3=="G"' \
  result/kraken2/bracken.all_levels.long.txt \
  > result/kraken2/bracken.G.from_all.long.txt

# 从long总表中提取种水平
awk -F '\t' 'NR==1 || $3=="S"' \
  result/kraken2/bracken.all_levels.long.txt \
  > result/kraken2/bracken.S.from_all.long.txt

# 从带注释总表中提取门水平relative表
awk -F '\t' 'NR==1 || $2=="P"' \
  result/kraken2/bracken.all_levels.relative.with_taxonomy.txt \
  > result/kraken2/bracken.P.from_all.relative.with_taxonomy.txt

# 从带注释总表中提取属水平relative表
awk -F '\t' 'NR==1 || $2=="G"' \
  result/kraken2/bracken.all_levels.relative.with_taxonomy.txt \
  > result/kraken2/bracken.G.from_all.relative.with_taxonomy.txt

# 从带注释总表中提取种水平relative表
awk -F '\t' 'NR==1 || $2=="S"' \
  result/kraken2/bracken.all_levels.relative.with_taxonomy.txt \
  > result/kraken2/bracken.S.from_all.relative.with_taxonomy.txt

echo "Kraken2 + Bracken analysis finished."
#bracken改进20260430end————————————
#bracken改进
#1228 cpu01
#执行bracken
    mkdir -p temp/bracken
    # 测序数据长度，通常为150，早期有100/75/50/25
    readLen=150
    # 20%样本中存在才保留
    prop=0.2
	#定义kraken库
	type=pluspf
	threshold=10

for tax in D P G S;do
  for i in $(tail -n+2 result/metadata.txt | cut -f1);do
    # 添加tax到输出文件名，避免覆盖
    bracken -d "${db}/kraken2/${type}/" \
      -i "temp/kraken2/${i}.report" \
      -r "${readLen}" -l "${tax}" -t $cpu \
      -o "temp/bracken/${i}_${tax}.brk" \
      -w "temp/bracken/${i}_${tax}.report"
    # 检查命令是否成功执行
    if [ $? -ne 0 ]; then
      echo "Error processing sample ${i} at tax level ${tax}" >&2
    fi
  done
done

    # 需要确认行数一致才能按以下方法合并      
    wc -l temp/bracken/*.report


#汇总结果
db=/home/tang/db
# 1. 生成 *.count 文件（无修改，完美）
input_dir="temp/bracken"
for tax in D P G S; do
  for i in $(tail -n+2 result/metadata.txt | cut -f1); do
    brk_file="${input_dir}/${i}_${tax}.brk"
    out_count="${input_dir}/${i}_${tax}.count"
    tail -n +2 "${brk_file}" | cut -f 6 | sed "1 s/^/${i}\n/" > "${out_count}"
    echo "已生成：${out_count}"
  done
done

# 2. 汇总合并（精简冗余，更稳定）
output_dir="result/kraken2"
prop=0.2

for tax in D P G S; do
  echo -e "\n===== 处理 $tax 级别 ====="

  # 精简：直接取第一个样本
  first_sample=$(tail -n+2 result/metadata.txt | cut -f1 | head -n1)
  brk_file="${input_dir}/${first_sample}_${tax}.brk"
  [ ! -f "$brk_file" ] && echo "跳过：$brk_file 不存在" && continue

  # 生成分类名表头
  tail -n+2 "$brk_file" | LC_ALL=C sort | cut -f1 | sed "1i Taxonomy" > "${input_dir}/0header_${tax}.count"
  [ ! -s "${input_dir}/0header_${tax}.count" ] && continue

  # 匹配所有样本count
  count_files=$(ls -1 "${input_dir}"/*_${tax}.count | sort)
  paste "${input_dir}/0header_${tax}.count" $count_files > "${output_dir}/bracken.${tax}.txt"

  # 统计
  csvtk -t stat "${output_dir}/bracken.${tax}.txt"

  # 过滤（你已成功的路径）
  Rscript "$db/EasyMicrobiome/script/filter_feature_table.R" \
    -i "${output_dir}/bracken.${tax}.txt" \
    -p "$prop" \
    -o "${output_dir}/bracken.${tax}.${prop}"

  # 清理临时文件
  rm -f "${input_dir}/0header_${tax}.count"
  echo "完成：${output_dir}/bracken.${tax}.${prop}"
done
#bracken改进end_______________________

#### 多样性和可视化
#alpha多样性计算：Berger Parker’s (BP), Simpson’s (Si), inverse Simpson’s (ISi), Shannon’s (Sh) # Fisher’s (F)依赖scipy.optimize包，默认未安装
    echo -e "SampleID\tBerger Parker\tSimpson\tinverse Simpson\tShannon" > result/kraken2/alpha.txt
    for i in `tail -n+2 result/metadata.txt|cut -f1`;do
        echo -e -n "$i\t" >> result/kraken2/alpha.txt
        for a in BP Si ISi Sh;do
            alpha_diversity.py -f temp/bracken/${i}.brk -a $a | cut -f 2 -d ':' | tr '\n' '\t' >> result/kraken2/alpha.txt
        done
        echo "" >> result/kraken2/alpha.txt
    done
    cat result/kraken2/alpha.txt
#beta多样性计算
    beta_diversity.py -i temp/bracken/*.brk --type bracken \
      > result/kraken2/beta.txt
    cat result/kraken2/beta.txt
#Krona图
    for i in `tail -n+2 result/metadata.txt|cut -f1`;do
        kreport2krona.py -r temp/bracken/${i}.report -o temp/bracken/${i}.krona --no-intermediate-ranks
        ktImportText temp/bracken/${i}.krona -o result/kraken2/krona.${i}.html
    done
#改进多样性和可视化

input_dir="temp/bracken"
output_dir="result/kraken2"
metadata_file="${output_dir%/*}/metadata.txt"

# 创建输出目录
mkdir -p "${output_dir}"

# 定义分类级别
tax_levels="D P G S"

# ====================== Alpha多样性计算（修复值提取逻辑） ======================
for tax in ${tax_levels}; do
    alpha_out="${output_dir}/alpha_${tax}.txt"
    echo -e "SampleID\tBerger Parker\tSimpson\tinverse Simpson\tShannon" > "${alpha_out}"

    # 读取样本ID
    sample_ids=$(tail -n+2 "${metadata_file}" | cut -f1 | sort)
    [ -z "${sample_ids}" ] && echo "ERROR: 未从metadata中提取到样本ID" >&2 && exit 1

    for i in ${sample_ids}; do
        brk_file="${input_dir}/${i}_${tax}.brk"
        [ ! -f "${brk_file}" ] && echo "WARNING: ${brk_file}不存在，跳过样本${i}(${tax}级)" >&2 && continue

        # 初始化一行的结果，先写入样本ID
        line_content="${i}"

        # 循环计算每个指数，修复值提取逻辑
        for a in BP Si ISi Sh; do
            # 先捕获原始输出，便于调试
            raw_output=$(alpha_diversity.py -f "${brk_file}" -a "${a}" 2>/dev/null)
            
            # 方法1：提取冒号后所有内容（处理空格/制表符）
            alpha_value=$(echo "${raw_output}" | awk -F ':' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}')
            
            # 方法2：如果方法1失败，直接提取数字（兜底方案）
            if [ -z "${alpha_value}" ] || [ "${alpha_value}" = "" ]; then
                alpha_value=$(echo "${raw_output}" | grep -o '[0-9]*\.[0-9]*\|[0-9]*' | head -n1)
            fi
            
            # 最终兜底为NA
            [ -z "${alpha_value}" ] && alpha_value="NA"

            # 拼接当前指数值（制表符分隔）
            line_content="${line_content}\t${alpha_value}"
        done

        # 将整行写入文件
        echo -e "${line_content}" >> "${alpha_out}"
    done

    echo -e "\n=== ${tax}级Alpha多样性计算完成 ==="
    # 只打印前5行预览，避免输出过长
    head -n5 "${alpha_out}"
    echo "..."
done

# ====================== Beta多样性计算（按分类级别） ======================
for tax in ${tax_levels}; do
    echo -e "\n=== 开始计算${tax}级Beta多样性 ==="
    brk_files=$(ls -1 "${input_dir}"/*_${tax}.brk 2>/dev/null | sort)
    [ -z "${brk_files}" ] && echo "WARNING: 未找到${input_dir}下的*_${tax}.brk文件，跳过${tax}级Beta多样性计算" >&2 || {
        beta_out="${output_dir}/beta_${tax}.txt"
        beta_diversity.py -i ${brk_files} --type bracken > "${beta_out}"
        echo "${tax}级Beta多样性计算完成，结果文件：${beta_out}"
        # 预览前5行
        head -n5 "${beta_out}"
        echo "..."
    }
done

# ====================== Krona图生成（按分类级别） ======================
for tax in ${tax_levels}; do
    echo -e "\n=== 开始生成${tax}级Krona图 ==="
    sample_ids=$(tail -n+2 "${metadata_file}" | cut -f1 | sort)
    for i in ${sample_ids}; do
        report_file="${input_dir}/${i}_${tax}.report"
        [ ! -f "${report_file}" ] && echo "WARNING: ${report_file}不存在，跳过样本${i}(${tax}级)的Krona图生成" >&2 && continue

        krona_temp="${input_dir}/${i}_${tax}.krona"
        kreport2krona.py -r "${report_file}" -o "${krona_temp}" --no-intermediate-ranks
        
        [ ! -s "${krona_temp}" ] && echo "WARNING: ${krona_temp}生成失败或为空，跳过样本${i}(${tax}级)" >&2 && continue

        krona_out="${output_dir}/krona.${i}_${tax}.html"
        ktImportText "${krona_temp}" -o "${krona_out}"
        [ $? -eq 0 ] && echo "✅ 样本${i}(${tax}级)：${krona_out}" || echo "❌ 样本${i}(${tax}级)生成失败" >&2

        rm -f "${krona_temp}"
    done
done


#改进多样性和可视化end————

#1228 cpu01end____________________________________________
#argoap
#argoap
#1224
wd=/home/tang/cssd
cpu=80
cd $wd
conda activate argoap
time args_oap stage_one -i $wd/temp/hr -o $wd/temp/argoap/output -f fastq -t $cpu  # 第一阶段，生成元数据等
time args_oap stage_two -i $wd/temp/argoap/output -t $cpu  # 第二阶段，生成归一化结果
#1224end——————————————



# 三、组装分析流程 Assemble-based

##  混合组装
20251228
    # 启动工作环境
    conda activate megahit
      # 小组ID: A1/A2/A3
    g=A1
g=A2
g=A3	
### MEGAHIT组装Assembly
#1972GB 需要k21需要1360GB内存；实际需要内存超过1.9tb
    # 删除旧文件夹，否则megahit无法运行
    # /bin/rm -rf temp/megahit
    # 组装，10~30m，32p18s8h, TB级数据需几天至几周，MEGAHIT v1.2.9
	#输入文件为temp/hr/*_1.fastq和*_2.fastq
	g=A1  
    time megahit -t $cpu \
        -1 `tail -n+2 result/metadata${g}.txt|cut -f1|sed 's/^/temp\/hr\//;s/$/_1.fastq/'|tr '\n' ','|sed 's/,$//'` \
        -2 `tail -n+2 result/metadata${g}.txt|cut -f1|sed 's/^/temp\/hr\//;s/$/_2.fastq/'|tr '\n' ','|sed 's/,$//'` \
        -o temp/megahit$g #--k-list 97,119 这个是可选，K值越大越快
		
    # 统计大小通常300M~5G，如18s100G10h1.8G
    # 如果contigs太多，可以按长度筛选，降低数据量，提高基因完整度，详见附录megahit
    seqkit stat temp/megahit$g/final.contigs.fa
    # 预览重叠群最前6行，前60列字符
    head -n6 temp/megahit$g/final.contigs.fa | cut -c1-60
	
	#可选删除低质量contig
	# megahit默认>200，可选 > 500 / 1000 bp，并统计前后变化；如此处筛选 > 500 bp，序列从15万变为3.5万条，总长度从7M下降到3M
    # mv temp/megahit/final.contigs.fa temp/megahit/raw.contigs.fa
    # seqkit seq -m 500 temp/megahit/raw.contigs.fa > temp/megahit/final.contigs.fa
    # seqkit stat temp/megahit/raw.contigs.fa
    # seqkit stat temp/megahit/final.contigs.fa
	#删除低质量contig end---------------------------
    # 备份重要结果
    mkdir -p result/megahit$g/
    ln -f temp/megahit$g/final.contigs.fa result/megahit$g/
    # 删除临时文件
    /bin/rm -rf temp/megahit$g/intermediate_contigs
	
20251228end_________________________

#===================================================================================================================================
20251229 09：30 1230 04：25  cpu01 rush 3 内存153GB
## 单样本组装与分箱
## (可选Opt)单样本分箱Single sample binning
tmux new -s cssdsb

    wd=/home/tang/cssd
	cpu=80
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
    mkdir -p seq temp result
	mkdir -p temp/megahit
    # p:threads线程数,job任务数,complete完整度x:contaminate污染率
p=80                     # 线程数
j=3                      # 并行任务数（rush的并发数）
c=50                     # 分箱提纯的完整性阈值（50%）
x=10                     # 分箱提纯的污染率阈值（10%）
#/home/tang/EasyMicrobiome/linux/rush
    
    conda activate metawrap    
    time tail -n+2 result/metadata2.txt|cut -f1|rush -j ${j} \
      "metawrap assembly -m 120 -t ${p} --megahit \
        -1 temp/hr/{}_1.fastq -2 temp/hr/{}_2.fastq \
        -o temp/megahit/{}"

#**分箱binning**
mkdir -p temp/binning
    time tail -n+2 result/metadata2.txt|cut -f1|rush -j ${j} \
      "metawrap binning \
        -o temp/binning/{} -t ${p} \
        -a temp/megahit/{}/final_assembly.fasta \
        --metabat2 --maxbin2 --concoct \
        temp/hr/{}_*.fastq > /dev/null 2>&1" 
20251230end_________________________
#**分箱提纯bin refinement**
20251231
mkdir -p temp/bin_refinement
    time tail -n+2 result/metadata2.txt|cut -f1|rush -j ${j} \
      "metawrap bin_refinement \
      -o temp/bin_refinement/{} -t ${p} \
      -A temp/binning/{}/metabat2_bins/ \
      -B temp/binning/{}/maxbin2_bins/ \
      -C temp/binning/{}/concoct_bins/ \
      -c ${c} -x ${x} "
    # 分别为1,2,2个
    tail -n+2 result/metadata2.txt|cut -f1|rush -j 1 \
      "tail -n+2 temp/bin_refinement/{}/metawrap_50_10_bins.stats|wc -l "
20251231enddone_________________________
#单样品分箱链接和重命名
pwd
wd=/home/tang/cssd
cd $wd

mkdir -p $wd/temp/drep_in
tail -n+2 result/metadata2.txt | cut -f1 | while read -r i; do
    # 检查样本目录是否存在
    bin_dir="$wd/temp/bin_refinementtest/${i}/metawrap_50_10_bins"
    if [ ! -d "$bin_dir" ]; then
        echo "警告：样本${i}的目录${bin_dir}不存在，跳过"
        continue
    fi
    
    # 遍历该样本下的所有bin文件
    for bin_file in "$bin_dir"/bin.*; do
        # 只处理实际存在的文件（避免匹配不到时的空值）
        if [ -f "$bin_file" ] || [ -L "$bin_file" ]; then
            # 获取bin文件名（如bin.1.fa）
            bin_name=$(basename "$bin_file")
            # 定义链接的临时路径（先创建原始名链接）
            temp_link="temp/drep_in/${bin_name}"
            # 定义最终重命名后的路径（添加Sg_${i}_前缀）
            final_name="temp/drep_in/Sg_${i}_${bin_name}"
            
            # 创建符号链接（-f强制覆盖已存在的链接）
            ln -sf "$bin_file" "$temp_link"
            # 用mv重命名（替代rename命令）
            mv -f "$temp_link" "$final_name"
            
            echo "成功处理：$bin_file -> $final_name"
        fi
    done
done
    # 删除空白中无效链接
    /bin/rm -f temp/drep_in/*\*
    # 统计混合和单样本来源数据，10个混，5个单；不同系统结果略有差异
    ls temp/drep_in/|cut -f 1 -d '_'|uniq -c
    # 统计混合批次/单样本来源
    ls temp/drep_in/|cut -f 2 -d '_'|cut -f 1 -d '.' |uniq -c
## 单样本组装与分箱enddone 测试完毕_______________________________________
20260227 cpu02 metadata2已经bin跑完，接着跑metadata1的bin再合并
tmux new -s cssdsb
    wd=/home/tang/cssd
	cpu=80
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
    mkdir -p seq temp result
	mkdir -p temp/megahit
    # p:threads线程数,job任务数,complete完整度x:contaminate污染率
p=80                     # 线程数
j=3                      # 并行任务数（rush的并发数）
c=50                     # 分箱提纯的完整性阈值（50%）
x=10                     # 分箱提纯的污染率阈值（10%）
#/home/tang/EasyMicrobiome/linux/rush
    conda activate metawrap   
#**分箱binning**
    time tail -n+2 result/metadata1.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j ${j} \
      "metawrap binning \
        -o temp/binningtest/{} -t ${p} \
        -a temp/megahit/{}/final_assembly.fasta \
        --metabat2 --maxbin2 --concoct \
        temp/hr/{}_*.fastq > /dev/null 2>&1" 

#**分箱提纯bin refinement**
    time tail -n+2 result/metadata1.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j ${j} \
      "metawrap bin_refinement \
      -o temp/bin_refinementtest/{} -t ${p} \
      -A temp/binningtest/{}/metabat2_bins/ \
      -B temp/binningtest/{}/maxbin2_bins/ \
      -C temp/binningtest/{}/concoct_bins/ \
      -c ${c} -x ${x} "
    # 分别为1,2,2个
    tail -n+2 result/metadata1.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j 1 \
      "tail -n+2 temp/bin_refinementtest/{}/metawrap_50_10_bins.stats|wc -l "
20260227 cpu02 end---------------------------
20260305 全部重新bin和refine

p=38                     # 线程数
j=2                      # 并行任务数（rush的并发数）
c=50                     # 分箱提纯的完整性阈值（50%）
x=10                     # 分箱提纯的污染率阈值（10%）
#/home/tang/EasyMicrobiome/linux/rush
    conda activate metawrap   
cd /home/tang/cssd
mkdir -p temp/binningtest
    time tail -n+2 result/metadata.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j ${j} \
      "metawrap binning \
        -o temp/binningtest/{} -t ${p} \
        -a result/megahit/{}_final_assembly.fasta \
        --metabat2 --maxbin2 --concoct \
        temp/hr/{}_*.fastq > /dev/null 2>&1" 
20260326
    time tail -n+2 result/metadata4.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j ${j} \
      "metawrap binning \
        -o temp/binningtest/{} -t ${p} \
        -a result/megahit/{}_final_assembly.fasta \
        --metabat2 --maxbin2 --concoct \
        temp/hr/{}_*.fastq > /dev/null 2>&1" 

#**分箱提纯bin refinement**
    mkdir -p temp/bin_refinementtest
	time tail -n+2 result/metadata.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j ${j} \
      "metawrap bin_refinement \
      -o temp/bin_refinementtest/{} -t ${p} \
      -A temp/binningtest/{}/metabat2_bins/ \
      -B temp/binningtest/{}/maxbin2_bins/ \
      -C temp/binningtest/{}/concoct_bins/ \
      -c ${c} -x ${x} "
    # 分别为1,2,2个
    tail -n+2 result/metadata1.txt|cut -f1|/home/tang/EasyMicrobiome/linux/rush -j 1 \
      "tail -n+2 temp/bin_refinementtest/{}/metawrap_50_10_bins.stats|wc -l "
20260330 done ，但是仅有SZ样本没有bin结果，重新处理	  
20260330
wd=/home/tang/cssd
cpu=76
soft=/home/tang/miniconda3
db=/home/tang/db
sample=SZ

mkdir -p $wd && cd $wd
mkdir -p seq temp result
mkdir -p temp/megahit temp/binning temp/bin_refinement

conda activate metawrap

########################################
# 1. 重新组装 assembly（仅 SZ）
########################################
rm -rf temp/megahit/${sample}

metawrap assembly -m 120 -t ${cpu} --megahit \
  -1 temp/hr/${sample}_1.fastq \
  -2 temp/hr/${sample}_2.fastq \
  -o temp/megahit/${sample}

########################################
# 2. 重新分箱 binning（仅 SZ）
########################################
rm -rf temp/binning/${sample}

metawrap binning \
  -o temp/binning/${sample} \
  -t ${cpu} \
  -a temp/megahit/${sample}/final_assembly.fasta \
  --metabat2 --maxbin2 --concoct \
  temp/hr/${sample}_1.fastq temp/hr/${sample}_2.fastq

########################################
# 3. 重新分箱提纯 bin refinement（仅 SZ）
########################################
rm -rf temp/bin_refinement/${sample}

metawrap bin_refinement \
  -o temp/bin_refinement/${sample} \
  -t ${cpu} \
  -A temp/binning/${sample}/metabat2_bins/ \
  -B temp/binning/${sample}/maxbin2_bins/ \
  -C temp/binning/${sample}/concoct_bins/ \
  -c 50 -x 10
20260330 done cpu07
=====================================================================================================================
!!!!!!组装和分箱过程！！！！！！！！
20251231 16：07 cpu01 rush 3 内存153GB
## 单样本组装与混合分箱
## (可选Opt)单样本分箱Single sample binning
tmux new -s cssdtest

    wd=/home/tang/cssd
	cpu=80
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
    mkdir -p seq temp result
	mkdir -p temp/megahit
    # p:threads线程数,job任务数,complete完整度x:contaminate污染率
p=80                     # 线程数
j=3                      # 并行任务数（rush的并发数）
c=50                     # 分箱提纯的完整性阈值（50%）
x=10                     # 分箱提纯的污染率阈值（10%）
    conda activate metawrap    
    time tail -n+2 result/metadata1.txt|cut -f1|rush -j ${j} \
      "metawrap assembly -m 120 -t ${p} --megahit \
        -1 temp/hr/{}_1.fastq -2 temp/hr/{}_2.fastq \
        -o temp/megahit/{}"
20260107 6：07 enddone——————————————————————————————————		

#！！！！！！！！！后续需要将所有单样本组装结果，先复制，再重命名，再合并！！！！！！！！！！！！！！！		

#将组装的结果复制——cpu07
for d in /home/tang/cssd/temp/megahit/*; do
    name=$(basename "$d")
    cp "$d/final_assembly.fasta" \
       "/home/tang/cssd/result/megahit/${name}_final_assembly.fasta"
done
#查看结果
ls -lh /home/tang/cssd/result/megahit

#检查contig序列命名是否重复
grep -h '^>' ~/cssd/result/megahit/*_final_assembly.fasta \
| sed 's/^>//' \
| awk '{print $1}' \
| sort \
| uniq -c \
| sort -nr \
| head
      2 k141_85890_length_1042_cov_4.0000
      2 k141_600481_length_1174_cov_3.0000
      2 k141_2795536_length_1308_cov_3.0000
      2 k141_2448219_length_1187_cov_3.0000
      2 k141_2447914_length_1017_cov_3.0000
      2 k141_1905760_length_1356_cov_4.0000
      1 k141_99999_length_3582_cov_5.0000
      1 k141_999996_length_1620_cov_4.0000
      1 k141_999995_length_1790_cov_6.0000
      1 k141_99998_length_1797_cov_4.0000
#以上结果表示重复

#给所有contig序列名称添加单独抬头
for f in ~/cssd/result/megahit/*_final_assembly.fasta; do
    prefix=$(basename "$f" _final_assembly.fasta)
    sed -i "s/^>/>${prefix}_/" "$f"
done

#将所有contig合并
cat ~/cssd/result/megahit/*_final_assembly.fasta \
> ~/cssd/result/megahit/all_samples_contigs.fasta

20260107 cpu01
mkdir -p temp/binning

time metawrap binning \
        -o temp/binning -t ${p} \
        -a result/megahit/all_samples_contigs.fasta \
        --metabat2 --maxbin2 --concoct \
        temp/hr/*.fastq > /dev/null 2>&1
20260113 cpu01 end——————————————————————————————————	


	
20260113 cpu01
    metawrap bin_refinement \
      -o temp/bin_refinement \
      -A temp/binning/metabat2_bins/ \
      -B temp/binning/maxbin2_bins/ \
      -c 70 -x 10 -t ${p}
20260113 cpu01 end  报错 了上个程序没跑完——————————————————————————————————

重跑
20260114 cpu01
time metawrap binning \
  -o temp/binning_l500 \
  -t ${p} \
  -l 500 \
  -a result/megahit/all_samples_contigs.fasta \
  --maxbin2 --concoct \
  temp/hr/*.fastq
20260114 cpu01end____________________________________________添加了关键参数 l 500，上次报错可能是应为contig默认设置的太短，导致maxbin和concoct没出结果，！！！！不对，默认是1000，设置成500更短了，导致结果更不对


还是报错，可能需要降低核心数量
20260120
p=64
time metawrap binning \
        -o temp/binning -t ${p} \
        -a result/megahit/all_samples_contigs.fasta \
        --metabat2 --maxbin2 --concoct \
        temp/hr/*.fastq
20260123end————————————————————————————
结果是metabat2跑完了，但是maxbin2出错
20260123
p=64
minlen=2500   # 2500/3000/5000 都可以试

time metawrap binning \
  -o temp/binning -t ${p} -m 256 \
  -a result/megahit/all_samples_contigs.fasta \
  -l ${minlen} \
  --metabat2 --maxbin2 --concoct \
  temp/hr/*.fastq
20260123
=====================================================================================================================
!!!!!!病毒组分析过程！！！！！！！！
#正式运行病毒的工具
# 1 genomad____________________________________
0107 10：15 cpu02   cssdvs2
conda activate genomad
time genomad end-to-end \
  ~/cssd/result/megahit/all_samples_contigs.fasta \
  ~/cssd/result/genomad \
  /home/tang/db/genomad_db \
  --threads 76
  
0107 10：15 cpu02 _____1231 06：25 cpu07end done——————————————————————————————————————————

#2 virsorter————————————————————————————————————————————————————————————————————————

0112  cpu02
 # 1. 启动代理（单独一个 terminal 或 tmux）-其实可以不用启动
 #/home/tang/clash/clash -f /home/tang/clash/1760165527084.yml

tmux new -s cssdvs2
conda activate vs2
time virsorter run \
  -w /home/tang/cssd/result/virsorter_out \
  -i /home/tang/cssd/result/megahit/all_samples_contigs.fasta \
  --min-length 1500 \
  --include-groups dsDNAphage,ssDNA,RNA,NCLDV \
  --db-dir /home/tang/db/virsorter/db \
  -j 76 \
  all
0112————enddone 跑了一个月

  # 1230 22：40 cpu05end_________________________________________
#仅保留
🧪 如果你只想重新筛选（以后会用到）

当你第一次 all 跑完后，可以只跑 classify（非常快）：

virsorter run classify \
  -w ~/cssd/result/virsorter2_all \
  --min-score 0.7 \
  --high-confidence-only
 #

#3.vibrant——————————————————————————————————————————————————————————————————
20260207 cpu04
tmux new -s cssdvib
conda activate vibrant
# time python3 /home/tang/mytools/VIBRANT/VIBRANT_run.py \  # 使用 python3 调用 VIBRANT 主程序
  # -i ~/cssd/result/megahit/all_samples_contigs.fasta \  # 输入文件：MEGAHIT 组装得到的所有样品 contigs（核酸序列）
  # -f nucl \                                           # 输入序列类型：核酸序列（nucl），而非蛋白序列
  # -folder ~/cssd/result/VIBRANT \                      # 输出目录：存放所有 VIBRANT 结果和临时文件
  # -t 76 \                                              # 并行运行的 scaffold 数量（一次同时分析 32 条 contig，每条占用 1 CPU）
  # -l 1500 \                                            # 最小 contig 长度阈值（bp），仅分析长度 ≥1500 bp 的序列
  # -o 4                                                 # 每条 contig 至少包含 4 个 ORFs 才参与病毒判定（默认值，平衡灵敏度和准确性）
time python3 /home/tang/mytools/VIBRANT/VIBRANT_run.py \
  -i ~/cssd/result/megahit/all_samples_contigs.fasta \
  -f nucl \
  -folder ~/cssd/result/VIBRANT \
  -t 76 \
  -l 1500 \
  -o 4
time python3 /home/tang/mytools/VIBRANT/VIBRANT_run.py \
  -i ~/cssd/result/megahit/all_samples_contigs.fasta \
  -f nucl \
  -folder ~/cssd/result/VIBRANT1 \
  -t 76 \
  -l 1500 \
  -o 4
报错，被killed可能是内存或者cpu过高,注意，消耗峰值内存1Tb以上
20260207 cpu07 end——————————————————————————————————

20260208 cpu07
time python3 /home/tang/mytools/VIBRANT/VIBRANT_run.py \
  -i ~/cssd/result/megahit/all_samples_contigs.fasta \
  -f nucl \
  -folder ~/cssd/result/VIBRANT2 \
  -t 76 \
  -l 1500 \
  -o 4 
20260208 cpu07 end  


# 4. cenote-taker3
20260208 cpu07
conda activate ct3_env
mkdir -p /home/tang/cssd/result/ct3_env

export CENOTE_DBS=/home/tang/db/ct3_dbs #添加数据库
time cenotetaker3 \
  -c /home/tang/cssd/result/megahit/all_samples_contigs.fasta \
  -r ct3_run \
  -p True \
  -wd /home/tang/cssd/result/ct3_env \
  -t 76 \
  --minimum_length_circular 1000 \
  --minimum_length_linear 1000

cenotetaker3 \
  -c /home/tang/cssd/result/megahit/all_samples_contigs.fasta \  # 输入contigs文件（必填）
  -r ct3_run \  # 运行名称（会作为输出子目录名，必填）
  -p True \  # 是否修剪原噬菌体区域（必填）
  -wd /home/tang/cssd/result/ct3_env \  # 工作目录（替代原-o参数，结果会存于此目录下的ct3_run子目录）
  -t 76 \  # 线程数
  --minimum_length_circular 1000 \  # 环形序列最小长度
  --minimum_length_linear 1000  # 线性序列最小长度
20260208 cpu07 end——————————————————————————————————
# 5. checkmV
20260212
conda activate checkv
# 核心运行命令：checkv end_to_end 输入faa文件 输出结果目录 数据库路径
checkv end_to_end /home/tang/cssd/result/megahit/all_samples_contigs.fasta /home/tang/cssd/result/checkv_result -d ~/db/checkv/checkv-db-v1.5 -t 76
20260212


=====================================================================================================================
20260427补充
wd=/home/tang/cssd
cpu=76
soft=/home/tang/miniconda3
db=/home/tang/db

mkdir -p $wd && cd $wd
mkdir -p seq temp result

############################################################
# 0. 目录与参数
############################################################
outdir=$wd/result/vOTU_formal_cdhit
tmpdir=$wd/temp/vOTU_formal_cdhit
mkdir -p $outdir $tmpdir

assembly=$wd/result/megahit/all_samples_contigs.fasta
reads_dir=$wd/temp/hr
checkv_db=$db/checkv/checkv-db-v1.5

# 已有病毒候选输出
genomad_virus=$wd/result/genomad/all_samples_contigs_summary/all_samples_contigs_virus.fna
genomad_provirus=$wd/result/genomad/all_samples_contigs_find_proviruses/all_samples_contigs_provirus.fna
vs2_virus=$wd/result/virsorter_out/final-viral-combined.fa
vibrant_virus=$wd/result/VIBRANT2/VIBRANT_all_samples_contigs/VIBRANT_phages_all_samples_contigs/all_samples_contigs.phages_combined.fna
ct3_virus=$wd/result/ct3_env/ct3_run/ct3_run_virus_sequences.fna

# 参数
map_cpu=$cpu

# 多工具支持阈值；更贴近参考论文可改成 1
MIN_SUPPORT=2

# 候选病毒序列进入聚类前的长度阈值
CANDIDATE_MINLEN=5000

# CD-HIT 近似 95% ANI + 85% AF
CDHIT_IDENTITY=0.95
CDHIT_COVERAGE_SHORTER=0.85

# CheckV 后导出高质量子集时使用
HQ_MINLEN=5000

############################################################
# 4. 统计各病毒工具 contig 支持数
############################################################
echo "[1/8] collecting contig ids from caller outputs ..."

caller_id_table=$tmpdir/caller_id_table.tsv

$soft/bin/conda run -n base --no-capture-output python - <<PY
from pathlib import Path

out = Path(r"$caller_id_table")
records = []

files = [
    ("geNomad_virus", r"$genomad_virus"),
    ("geNomad_provirus", r"$genomad_provirus"),
    ("VirSorter2", r"$vs2_virus"),
    ("VIBRANT", r"$vibrant_virus"),
    ("CT3", r"$ct3_virus"),
]

def fasta_ids(fp):
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]

seen = set()
for caller, fp in files:
    for sid in fasta_ids(fp):
        key = (caller, sid)
        if key not in seen:
            records.append((caller, sid))
            seen.add(key)

with open(out, "w") as fo:
    for caller, sid in records:
        fo.write(f"{caller}\t{sid}\n")
PY

support_table=$outdir/contig_tool_support.tsv

$soft/bin/conda run -n base --no-capture-output python - <<PY
from collections import defaultdict

inp = r"$caller_id_table"
out = r"$support_table"

d = defaultdict(set)
with open(inp) as f:
    for line in f:
        caller, sid = line.rstrip().split("\t")
        d[sid].add(caller)

rows = []
for sid, callers in d.items():
    rows.append((sid, len(callers), ",".join(sorted(callers))))
rows.sort(key=lambda x: (-x[1], x[0]))

with open(out, "w") as fo:
    fo.write("contig_id\ttool_support_n\ttools\n")
    for sid, n, tools in rows:
        fo.write(f"{sid}\t{n}\t{tools}\n")

with open(r"$tmpdir/ids_support_ge${MIN_SUPPORT}.txt", "w") as fo:
    for sid, n, tools in rows:
        if n >= $MIN_SUPPORT:
            fo.write(sid + "\n")
PY

############################################################
# 5. 从 assembly 提取候选病毒 contigs，先限制为 >=5000 bp
############################################################
echo "[2/8] extracting candidate viral contigs from assembly ..."

$soft/bin/conda run -n base --no-capture-output python - <<PY
from pathlib import Path

assembly = Path(r"$assembly")
keep_ids = set(x.strip() for x in open(r"$tmpdir/ids_support_ge${MIN_SUPPORT}.txt") if x.strip())
outfa = Path(r"$outdir/viral_candidates_ge${MIN_SUPPORT}_tools_min${CANDIDATE_MINLEN}.fa")
stats = Path(r"$outdir/viral_candidates_ge${MIN_SUPPORT}_tools_min${CANDIDATE_MINLEN}.stats.tsv")
minlen = $CANDIDATE_MINLEN

def fasta_iter(fp):
    sid = None
    seq = []
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(seq)
                sid = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if sid is not None:
            yield sid, "".join(seq)

n = 0
total_bp = 0
with open(outfa, "w") as fo:
    for sid, seq in fasta_iter(assembly):
        if sid in keep_ids and len(seq) >= minlen:
            fo.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i+80] + "\n")
            n += 1
            total_bp += len(seq)

with open(stats, "w") as so:
    so.write("file\tnseq\ttotal_bp\n")
    so.write(f"{outfa.name}\t{n}\t{total_bp}\n")
PY

candidate_fa=$outdir/viral_candidates_ge${MIN_SUPPORT}_tools_min${CANDIDATE_MINLEN}.fa

############################################################
# 6. 先用 CD-HIT 聚类，获得预 vOTU 代表序列
############################################################
echo "[3/8] clustering candidates with CD-HIT before CheckV ..."

cdhit_rep=$outdir/vOTU_rep_preCheckV.fa
cdhit_clstr=$outdir/vOTU_rep_preCheckV.fa.clstr

$soft/bin/conda run -n megahit --no-capture-output cd-hit-est \
  -i "$candidate_fa" \
  -o "$cdhit_rep" \
  -c "$CDHIT_IDENTITY" \
  -aS "$CDHIT_COVERAGE_SHORTER" \
  -g 1 \
  -G 0 \
  -T "$cpu" \
  -M 0 \
  -d 0

# 解析 .clstr -> vOTU 成员表
$soft/bin/conda run -n base --no-capture-output python - <<PY
from pathlib import Path
import re

clstr = Path(r"$cdhit_clstr")
out_members = Path(r"$outdir/vOTU_members.tsv")
out_summary = Path(r"$outdir/vOTU_summary.tsv")
out_map = Path(r"$outdir/vOTU_rep_map.tsv")

clusters = []
current = None
rep_pat = re.compile(r'>\s*([^\.\s]+(?:\.[^\.\s]+)*)')

with open(clstr) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">Cluster"):
            if current is not None:
                clusters.append(current)
            current = {"members": [], "rep": None}
        else:
            m = rep_pat.search(line)
            if not m:
                continue
            sid = m.group(1)
            current["members"].append(sid)
            if line.endswith("*"):
                current["rep"] = sid

if current is not None:
    clusters.append(current)

with open(out_members, "w") as fo:
    fo.write("vOTU_id\trepresentative_id\tmember_id\tn_members\n")
    for i, clu in enumerate(clusters, start=1):
        vid = f"vOTU_{i:05d}"
        rep = clu["rep"] if clu["rep"] is not None else clu["members"][0]
        n = len(clu["members"])
        for m in clu["members"]:
            fo.write(f"{vid}\t{rep}\t{m}\t{n}\n")

with open(out_summary, "w") as fo:
    fo.write("vOTU_id\trepresentative_id\tn_members\n")
    for i, clu in enumerate(clusters, start=1):
        vid = f"vOTU_{i:05d}"
        rep = clu["rep"] if clu["rep"] is not None else clu["members"][0]
        n = len(clu["members"])
        fo.write(f"{vid}\t{rep}\t{n}\n")

with open(out_map, "w") as fo:
    fo.write("vOTU_id\trepresentative_id\n")
    for i, clu in enumerate(clusters, start=1):
        vid = f"vOTU_{i:05d}"
        rep = clu["rep"] if clu["rep"] is not None else clu["members"][0]
        fo.write(f"{vid}\t{rep}\n")
PY

############################################################
# 7. 对聚类后的 vOTU 代表序列跑 CheckV
############################################################
echo "[4/8] running CheckV on clustered vOTU representatives ..."

rm -rf $outdir/checkv_on_vOTU
$soft/bin/conda run -n checkv --no-capture-output checkv end_to_end \
  "$cdhit_rep" \
  "$outdir/checkv_on_vOTU" \
  -d "$checkv_db" \
  -t "$cpu"

$soft/bin/conda run -n base --no-capture-output python - <<PY
import csv
from pathlib import Path

qsum = Path(r"$outdir/checkv_on_vOTU/quality_summary.tsv")
repmap = Path(r"$outdir/vOTU_rep_map.tsv")
out = Path(r"$outdir/vOTU_checkv_summary.tsv")

rep_to_votu = {}
with open(repmap) as f:
    next(f)
    for line in f:
        votu_id, rep = line.rstrip().split("\t")
        rep_to_votu[rep] = votu_id

with open(qsum) as f, open(out, "w") as fo:
    reader = csv.DictReader(f, delimiter="\t")
    fo.write("vOTU_id\trepresentative_id\tcontig_length\tquality\tcompleteness\tcontamination\n")
    for r in reader:
        rep = r["contig_id"]
        votu_id = rep_to_votu.get(rep, "NA")
        quality = r.get("checkv_quality", r.get("miuvig_quality", "NA"))
        completeness = r.get("completeness", "NA")
        contamination = r.get("contamination", "NA")
        contig_length = r.get("contig_length", "NA")
        fo.write(f"{votu_id}\t{rep}\t{contig_length}\t{quality}\t{completeness}\t{contamination}\n")
PY

############################################################
# 8. 可选：导出中高质量 vOTU 子集
############################################################
echo "[5/8] preparing optional HQ vOTU subset ..."

$soft/bin/conda run -n base --no-capture-output python - <<PY
from pathlib import Path
import csv

qsum = Path(r"$outdir/checkv_on_vOTU/quality_summary.tsv")
repfa = Path(r"$cdhit_rep")
keep_ids = Path(r"$tmpdir/votu_hq_ids.txt")
outfa = Path(r"$outdir/vOTU_rep_HQ.fa")
outsum = Path(r"$outdir/vOTU_checkv_HQ.tsv")
minlen = $HQ_MINLEN

keep = set()
rows = []

with open(qsum) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        rep = r["contig_id"]
        clen = int(r["contig_length"])
        quality = r.get("checkv_quality", r.get("miuvig_quality", ""))
        if quality in ("Complete", "High-quality", "Medium-quality") and clen >= minlen:
            keep.add(rep)
            rows.append((rep, clen, quality, r.get("completeness", "NA"), r.get("contamination", "NA")))

with open(keep_ids, "w") as fo:
    for rep in sorted(keep):
        fo.write(rep + "\n")

with open(outsum, "w") as fo:
    fo.write("representative_id\tcontig_length\tquality\tcompleteness\tcontamination\n")
    for x in rows:
        fo.write("\t".join(map(str, x)) + "\n")

def fasta_iter(fp):
    sid = None
    seq = []
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(seq)
                sid = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if sid is not None:
            yield sid, "".join(seq)

with open(outfa, "w") as fo:
    for sid, seq in fasta_iter(repfa):
        if sid in keep:
            fo.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i+80] + "\n")
PY

############################################################
# 9. 用全部聚类后的 vOTU 代表序列做丰度矩阵
############################################################
echo "[6/8] quantifying vOTU abundance with coverm ..."

abund_dir=$outdir/abundance_per_sample
mkdir -p $abund_dir

# metadata.txt 第一列为样本名；使用 rush 并行
# 如需调大并行数，可把 -j 2 改大；单任务线程数由 map_cpu 控制

tail -n+2 $wd/result/metadata.txt | cut -f1 | /home/tang/db/EasyMicrobiome/linux/rush -j 6 \
  "r1=$reads_dir/{}_1.fastq; \
   r2=$reads_dir/{}_2.fastq; \
   [[ -s \"\$r1\" ]] || r1=$reads_dir/{}_1.fastq.gz; \
   [[ -s \"\$r2\" ]] || r2=$reads_dir/{}_2.fastq.gz; \
   if [[ ! -s \"\$r1\" || ! -s \"\$r2\" ]]; then \
     echo 'WARN: missing reads for {}' > $abund_dir/{}.log; \
     exit 0; \
   fi; \
   $soft/bin/conda run -n coverm --no-capture-output coverm contig \
     --coupled \"\$r1\" \"\$r2\" \
     --reference $cdhit_rep \
     --methods tpm mean covered_fraction \
     --min-read-percent-identity 95 \
     --min-covered-fraction 0.75 \
     -t 10 \
     > $abund_dir/{}.tsv 2> $abund_dir/{}.log"


$soft/bin/conda run -n base --no-capture-output python - <<PY
import os, glob, re
import pandas as pd

abund_dir = r"$abund_dir"
files = sorted(glob.glob(os.path.join(abund_dir, "*.tsv")))
if not files:
    raise SystemExit("No abundance files found.")

def norm(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())

def find_col(cols, target):
    nt = norm(target)
    for c in cols:
        if nt == norm(c):
            return c
    for c in cols:
        if nt in norm(c):
            return c
    return None

tpm_list = []
mean_list = []
cov_list = []

for f in files:
    sample = os.path.basename(f).rsplit(".", 1)[0]
    df = pd.read_csv(f, sep="\t")
    idcol = df.columns[0]

    tpm_col = find_col(df.columns, "tpm")
    mean_col = find_col(df.columns, "mean")
    cov_col = find_col(df.columns, "covered_fraction")

    if tpm_col is None:
        raise SystemExit(f"TPM column missing in {f}")

    tpm_list.append(df.set_index(idcol)[tpm_col].rename(sample))

    if mean_col is not None:
        mean_list.append(df.set_index(idcol)[mean_col].rename(sample))

    if cov_col is not None:
        cov_list.append(df.set_index(idcol)[cov_col].rename(sample))

tpm = pd.concat(tpm_list, axis=1).fillna(0)
tpm.index.name = "representative_id"
tpm.to_csv(r"$outdir/vOTU_abundance_TPM.tsv", sep="\t")

if mean_list:
    mean = pd.concat(mean_list, axis=1).fillna(0)
    mean.index.name = "representative_id"
    mean.to_csv(r"$outdir/vOTU_abundance_Mean.tsv", sep="\t")

if cov_list:
    cov = pd.concat(cov_list, axis=1).fillna(0)
    cov.index.name = "representative_id"
    cov.to_csv(r"$outdir/vOTU_abundance_CoveredFraction.tsv", sep="\t")
PY


############################################################
# 10. 病毒携带 ARGs 识别（SARG + DIAMOND BLASTX）
############################################################
echo "[7/8] identifying viral ARGs with DIAMOND BLASTX ..."

mkdir -p $tmpdir/SARG $outdir/SARG

sarg_fasta=$db/SARG/4.SARG_v3.2_20220917_Short_subdatabase.fasta
sarg_db=$db/SARG/SARG
virus_arg_query=$outdir/vOTU_rep_HQ.fa

# BLASTX：核酸序列对蛋白数据库
conda activate eggnog
$soft/bin/conda run -n eggnog --no-capture-output diamond blastx \
  --db $sarg_db \
  --query $virus_arg_query \
  --threads 40 \
  -e 1e-7 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
  --max-target-seqs 10 \
  --quiet \
  --out $tmpdir/SARG/virus_ARG_diamond.tsv

# 过滤：identity >= 40%, hit length >= 75 aa, E-value <= 1e-7  未执行
$soft/bin/conda run -n base --no-capture-output Rscript - <<RS
suppressPackageStartupMessages(library(tidyverse))

inp  <- "$tmpdir/SARG/virus_ARG_diamond.tsv"
out1 <- "$outdir/SARG/virus_ARG_diamond.filtered.tsv"
out2 <- "$outdir/SARG/virus_ARG_diamond.filtered.best.tsv"
out3 <- "$outdir/SARG/virus_ARG_subject_count.tsv"
out4 <- "$outdir/SARG/vOTU_ARG_burden.tsv"

cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
)

if (!file.exists(inp) || file.info(inp)$size == 0) {
  write_tsv(tibble(
    qseqid = character(),
    sseqid = character(),
    pident = double(),
    length = integer(),
    evalue = double(),
    bitscore = double(),
    qlen = integer(),
    slen = integer(),
    coverage_query = double(),
    coverage_subject = double()
  ), out1)

  write_tsv(tibble(
    qseqid = character(),
    sseqid = character(),
    pident = double(),
    length = integer(),
    evalue = double(),
    bitscore = double(),
    qlen = integer(),
    slen = integer(),
    coverage_query = double(),
    coverage_subject = double()
  ), out2)

  write_tsv(tibble(
    sseqid = character(),
    n_vOTU = integer()
  ), out3)

  write_tsv(tibble(
    representative_id = character(),
    has_ARG = integer(),
    n_ARG_hit = integer()
  ), out4)

  quit(save = "no")
}

df <- read_tsv(
  inp,
  col_names = cols,
  show_col_types = FALSE,
  progress = FALSE
)

filtered <- df %>%
  mutate(
    pident   = as.numeric(pident),
    length   = as.integer(length),
    evalue   = as.numeric(evalue),
    bitscore = as.numeric(bitscore),
    qlen     = as.integer(qlen),
    slen     = as.integer(slen)
  ) %>%
  filter(
    pident >= 40,
    length >= 75,
    evalue <= 1e-7
  ) %>%
  mutate(
    coverage_query   = if_else(qlen > 0, round(length / qlen * 100, 3), 0),
    coverage_subject = if_else(slen > 0, round(length / slen * 100, 3), 0)
  ) %>%
  select(
    qseqid, sseqid, pident, length, evalue, bitscore,
    qlen, slen, coverage_query, coverage_subject
  )

write_tsv(filtered, out1)

best <- filtered %>%
  arrange(qseqid, desc(bitscore), desc(pident), desc(length), evalue) %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup()

write_tsv(best, out2)

subject_count <- best %>%
  count(sseqid, name = "n_vOTU")
write_tsv(subject_count, out3)

burden <- best %>%
  count(qseqid, name = "n_ARG_hit") %>%
  transmute(
    representative_id = qseqid,
    has_ARG = 1L,
    n_ARG_hit = n_ARG_hit
  )
write_tsv(burden, out4)
RS


############################################################
# 10. 输出
############################################################
echo "[7/8] done"
echo "Key outputs:"
echo "  $outdir/contig_tool_support.tsv"
echo "  $outdir/viral_candidates_ge${MIN_SUPPORT}_tools_min${CANDIDATE_MINLEN}.fa"
echo "  $outdir/vOTU_rep_preCheckV.fa"
echo "  $outdir/vOTU_rep_preCheckV.fa.clstr"
echo "  $outdir/vOTU_members.tsv"
echo "  $outdir/vOTU_summary.tsv"
echo "  $outdir/checkv_on_vOTU/"
echo "  $outdir/vOTU_checkv_summary.tsv"
echo "  $outdir/vOTU_rep_HQ.fa"
echo "  $outdir/vOTU_checkv_HQ.tsv"
echo "  $outdir/vOTU_abundance_TPM.tsv"
echo "  $outdir/vOTU_abundance_Mean.tsv"
echo "  $outdir/vOTU_abundance_CoveredFraction.tsv"

echo "  $outdir/vOTU_rep_HQ.fa"
这个是高质量病毒数据，可以用来跑ARGs注释
参考过滤参数BLASTX of DIAMOND v2.0.14.15254 was used to align the sequences, and a sequence was classified as an ARG fragment if it met the following criteria: identity ≥ 40%, hit length ≥ 75%, and E-value ≤ 10-7.




20260426 自选过滤方案  结果与chhit结果差不多
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

############################################################
# 0. 基本配置
############################################################
WD="/home/tang/cssd"
OUTDIR="$WD/result/vOTU_formal"
TMPDIR="$WD/temp/vOTU_formal"
mkdir -p "$OUTDIR" "$TMPDIR"

ASSEMBLY="$WD/result/megahit/all_samples_contigs.fasta"
READS_DIR="$WD/temp/hr"
CHECKV_DB="$HOME/db/checkv/checkv-db-v1.5"

# 已确认的病毒候选输出
GENOMAD_VIRUS="$WD/result/genomad/all_samples_contigs_summary/all_samples_contigs_virus.fna"
GENOMAD_PROVIRUS="$WD/result/genomad/all_samples_contigs_find_proviruses/all_samples_contigs_provirus.fna"
VS2_VIRUS="$WD/result/virsorter_out/final-viral-combined.fa"
VIBRANT_VIRUS="$WD/result/VIBRANT2/VIBRANT_all_samples_contigs/VIBRANT_phages_all_samples_contigs/all_samples_contigs.phages_combined.fna"
CT3_VIRUS="$WD/result/ct3_env/ct3_run/ct3_run_virus_sequences.fna"

# 参数
THREADS=76
MAP_THREADS=20
MIN_SUPPORT=2
PRECHECK_MINLEN=1500
POSTCHECK_MINLEN=5000
VOTU_ANI=95
VOTU_AF=85

############################################################
# 1. conda 环境配置
############################################################
MINICONDA_HOME="/home/tang/miniconda3"
source "$MINICONDA_HOME/etc/profile.d/conda.sh"

PY_ENV="base"
CHECKV_ENV="checkv"
BLAST_ENV="eggnog"
COVERM_ENV="coverm"

run_env () {
  local env_name="$1"
  shift
  conda run -n "$env_name" --no-capture-output "$@"
}

############################################################
# 2. 检查输入文件
############################################################
for f in \
  "$ASSEMBLY" \
  "$GENOMAD_VIRUS" \
  "$GENOMAD_PROVIRUS" \
  "$VS2_VIRUS" \
  "$VIBRANT_VIRUS" \
  "$CT3_VIRUS"
do
  [[ -s "$f" ]] || { echo "ERROR: missing file: $f"; exit 1; }
done

[[ -d "$READS_DIR" ]] || { echo "ERROR: missing dir: $READS_DIR"; exit 1; }
[[ -d "$CHECKV_DB" ]] || { echo "ERROR: missing CheckV DB: $CHECKV_DB"; exit 1; }

############################################################
# 3. 检查关键命令
############################################################
run_env "$PY_ENV" python -V >/dev/null
run_env "$CHECKV_ENV" checkv -h >/dev/null 2>&1 || { echo "ERROR: checkv not found in env $CHECKV_ENV"; exit 1; }
run_env "$BLAST_ENV" blastn -version >/dev/null 2>&1 || { echo "ERROR: blastn not found in env $BLAST_ENV"; exit 1; }
run_env "$BLAST_ENV" makeblastdb -version >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found in env $BLAST_ENV"; exit 1; }
run_env "$COVERM_ENV" coverm --help >/dev/null 2>&1 || { echo "ERROR: coverm not found in env $COVERM_ENV"; exit 1; }

############################################################
# 4. 统计每条 contig 被几个工具支持
############################################################
echo "[1/8] collecting contig ids from caller outputs ..."

CALLER_ID_TABLE="$TMPDIR/caller_id_table.tsv"

run_env "$PY_ENV" python - <<PY
from pathlib import Path

out = Path(r"$CALLER_ID_TABLE")
records = []

files = [
    ("geNomad_virus", r"$GENOMAD_VIRUS"),
    ("geNomad_provirus", r"$GENOMAD_PROVIRUS"),
    ("VirSorter2", r"$VS2_VIRUS"),
    ("VIBRANT", r"$VIBRANT_VIRUS"),
    ("CT3", r"$CT3_VIRUS"),
]

def fasta_ids(fp):
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]

seen = set()
for caller, fp in files:
    for sid in fasta_ids(fp):
        key = (caller, sid)
        if key not in seen:
            records.append((caller, sid))
            seen.add(key)

with open(out, "w") as fo:
    for caller, sid in records:
        fo.write(f"{caller}\t{sid}\n")
PY

SUPPORT_TABLE="$OUTDIR/contig_tool_support.tsv"

run_env "$PY_ENV" python - <<PY
from collections import defaultdict

inp = r"$CALLER_ID_TABLE"
out = r"$SUPPORT_TABLE"

d = defaultdict(set)
with open(inp) as f:
    for line in f:
        caller, sid = line.rstrip().split("\t")
        d[sid].add(caller)

rows = []
for sid, callers in d.items():
    rows.append((sid, len(callers), ",".join(sorted(callers))))
rows.sort(key=lambda x: (-x[1], x[0]))

with open(out, "w") as fo:
    fo.write("contig_id\ttool_support_n\ttools\n")
    for sid, n, tools in rows:
        fo.write(f"{sid}\t{n}\t{tools}\n")

with open(r"$TMPDIR/ids_support_ge${MIN_SUPPORT}.txt", "w") as fo:
    for sid, n, tools in rows:
        if n >= $MIN_SUPPORT:
            fo.write(sid + "\n")
PY

############################################################
# 5. 从总 assembly 提取多工具支持的 contig，并做长度过滤
############################################################
echo "[2/8] extracting multi-caller supported contigs from assembly ..."

run_env "$PY_ENV" python - <<PY
from pathlib import Path

assembly = Path(r"$ASSEMBLY")
keep_ids = set(x.strip() for x in open(r"$TMPDIR/ids_support_ge${MIN_SUPPORT}.txt") if x.strip())
outfa = Path(r"$OUTDIR/merged_candidates.precheck.fa")
stats = Path(r"$OUTDIR/merged_candidates.precheck.stats.tsv")
minlen = $PRECHECK_MINLEN

def fasta_iter(fp):
    sid = None
    seq = []
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(seq)
                sid = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if sid is not None:
            yield sid, "".join(seq)

n = 0
total_bp = 0
with open(outfa, "w") as fo:
    for sid, seq in fasta_iter(assembly):
        if sid in keep_ids and len(seq) >= minlen:
            fo.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i+80] + "\n")
            n += 1
            total_bp += len(seq)

with open(stats, "w") as so:
    so.write("file\tnseq\ttotal_bp\n")
    so.write(f"{outfa.name}\t{n}\t{total_bp}\n")
PY

############################################################
# 6. 对合并候选重跑 CheckV
############################################################
echo "[3/8] running CheckV ..."

rm -rf "$OUTDIR/checkv_final"
run_env "$CHECKV_ENV" checkv end_to_end \
  "$OUTDIR/merged_candidates.precheck.fa" \
  "$OUTDIR/checkv_final" \
  -d "$CHECKV_DB" \
  -t "$THREADS"

############################################################
# 7. 筛选 Complete / High-quality / Medium-quality
############################################################
echo "[4/8] filtering CheckV results ..."

run_env "$PY_ENV" python - <<PY
from pathlib import Path
import csv

qsum = Path(r"$OUTDIR/checkv_final/quality_summary.tsv")
virus_fa = Path(r"$OUTDIR/checkv_final/viruses.fna")
provirus_fa = Path(r"$OUTDIR/checkv_final/proviruses.fna")

kept_tsv = Path(r"$OUTDIR/checkv_kept.tsv")
kept_ids = Path(r"$TMPDIR/checkv_kept_ids.txt")
merged_fa = Path(r"$TMPDIR/checkv_all.fa")
outfa = Path(r"$OUTDIR/final_viral_hq.fa")
stats = Path(r"$OUTDIR/final_viral_hq.stats.tsv")

minlen = $POSTCHECK_MINLEN

def fasta_iter(fp):
    sid = None
    seq = []
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(seq)
                sid = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if sid is not None:
            yield sid, "".join(seq)

with open(merged_fa, "w") as fo:
    for fp in [virus_fa, provirus_fa]:
        if fp.exists() and fp.stat().st_size > 0:
            with open(fp) as fi:
                for line in fi:
                    fo.write(line)

keep = {}
with open(qsum) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        cid = r["contig_id"]
        clen = int(r["contig_length"])
        quality = r.get("checkv_quality", r.get("miuvig_quality", ""))
        comp = r.get("completeness", "NA")
        if quality in ("Complete", "High-quality", "Medium-quality") and clen >= minlen:
            keep[cid] = (clen, quality, comp)

with open(kept_tsv, "w") as fo:
    fo.write("contig_id\tcontig_length\tquality\tcompleteness\n")
    for cid, (clen, q, c) in sorted(keep.items()):
        fo.write(f"{cid}\t{clen}\t{q}\t{c}\n")

with open(kept_ids, "w") as fo:
    for cid in sorted(keep):
        fo.write(cid + "\n")

n = 0
total_bp = 0
with open(outfa, "w") as fo:
    for sid, seq in fasta_iter(merged_fa):
        if sid in keep:
            fo.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i+80] + "\n")
            n += 1
            total_bp += len(seq)

with open(stats, "w") as so:
    so.write("file\tnseq\ttotal_bp\n")
    so.write(f"{outfa.name}\t{n}\t{total_bp}\n")
PY

############################################################
# 8. all-vs-all blastn，按 95/85 聚成 vOTU
############################################################
echo "[5/8] all-vs-all blastn clustering ..."

run_env "$BLAST_ENV" makeblastdb \
  -in "$OUTDIR/final_viral_hq.fa" \
  -dbtype nucl \
  -out "$TMPDIR/final_viral_hq_db" >/dev/null

run_env "$BLAST_ENV" blastn \
  -task megablast \
  -query "$OUTDIR/final_viral_hq.fa" \
  -db "$TMPDIR/final_viral_hq_db" \
  -out "$TMPDIR/all_vs_all.blast.tsv" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  -evalue 1e-10 \
  -max_target_seqs 100000 \
  -num_threads "$THREADS"

run_env "$PY_ENV" python - <<PY
import csv
from collections import defaultdict

fasta = r"$OUTDIR/final_viral_hq.fa"
blast = r"$TMPDIR/all_vs_all.blast.tsv"
checkv_kept = r"$OUTDIR/checkv_kept.tsv"
ani_cutoff = float("$VOTU_ANI")
af_cutoff = float("$VOTU_AF")

lengths = {}
with open(fasta) as fh:
    sid = None
    seq = []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if sid is not None:
                lengths[sid] = sum(len(x) for x in seq)
            sid = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if sid is not None:
        lengths[sid] = sum(len(x) for x in seq)

meta = {}
with open(checkv_kept) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        cid = r["contig_id"]
        comp = r["completeness"]
        meta[cid] = {
            "quality": r["quality"],
            "completeness": float(comp) if comp not in ("NA", "", None) else -1.0
        }

def qrank(q):
    q = q.lower()
    if "complete" in q:
        return 3
    if "high-quality" in q:
        return 2
    if "medium-quality" in q:
        return 1
    return 0

def merge_intervals(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s + 1 for s, e in merged)

pair_hits = defaultdict(list)

with open(blast) as fh:
    for line in fh:
        q, s, pident, aln, mm, gap, qs, qe, ss, se, e, bits, qlen, slen = line.rstrip().split("\t")
        if q == s:
            continue
        pident = float(pident)
        aln = int(aln)
        qs, qe = sorted((int(qs), int(qe)))
        ss, se = sorted((int(ss), int(se)))
        key = tuple(sorted((q, s)))
        pair_hits[key].append({
            "q": q,
            "s": s,
            "pident": pident,
            "aln": aln,
            "qint": (qs, qe),
            "sint": (ss, se)
        })

pair_stats = {}
for (a, b), hits in pair_hits.items():
    qints = []
    sints = []
    ident_weighted = 0.0
    aln_total = 0

    for h in hits:
        if h["q"] == a and h["s"] == b:
            qints.append(h["qint"])
            sints.append(h["sint"])
        else:
            qints.append(h["sint"])
            sints.append(h["qint"])
        ident_weighted += h["pident"] * h["aln"]
        aln_total += h["aln"]

    qcov = merge_intervals(qints)
    scov = merge_intervals(sints)
    af = min(qcov / lengths[a], scov / lengths[b]) * 100.0
    ani = ident_weighted / aln_total if aln_total > 0 else 0.0

    pair_stats[(a, b)] = {
        "ani": ani,
        "af": af,
        "aln_total": aln_total
    }

parent = {k: k for k in lengths}
rank = {k: 0 for k in lengths}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        parent[ra] = rb
    elif rank[ra] > rank[rb]:
        parent[rb] = ra
    else:
        parent[rb] = ra
        rank[ra] += 1

for (a, b), st in pair_stats.items():
    if st["ani"] >= ani_cutoff and st["af"] >= af_cutoff:
        union(a, b)

clusters = defaultdict(list)
for sid in lengths:
    clusters[find(sid)].append(sid)

def rep_score(cid):
    q = meta.get(cid, {}).get("quality", "")
    comp = meta.get(cid, {}).get("completeness", -1.0)
    ln = lengths.get(cid, 0)
    return (qrank(q), comp, ln, cid)

cluster_items = []
for root, members in clusters.items():
    members = sorted(members)
    rep = sorted(members, key=rep_score, reverse=True)[0]
    cluster_items.append((rep, members))

cluster_items = sorted(cluster_items, key=lambda x: (len(x[1]), lengths.get(x[0], 0)), reverse=True)

with open(r"$OUTDIR/vOTU_edges.tsv", "w") as fo:
    fo.write("contig_a\tcontig_b\tani\taf_shorter\taln_total\n")
    for (a, b), st in sorted(pair_stats.items()):
        if st["ani"] >= ani_cutoff and st["af"] >= af_cutoff:
            fo.write(f"{a}\t{b}\t{st['ani']:.3f}\t{st['af']:.3f}\t{st['aln_total']}\n")

with open(r"$OUTDIR/vOTU_members.tsv", "w") as fo:
    fo.write("vOTU_id\trepresentative_id\tmember_id\tn_members\n")
    for i, (rep, members) in enumerate(cluster_items, start=1):
        vid = f"vOTU_{i:05d}"
        for m in members:
            fo.write(f"{vid}\t{rep}\t{m}\t{len(members)}\n")

with open(r"$OUTDIR/vOTU_summary.tsv", "w") as fo:
    fo.write("vOTU_id\trepresentative_id\tn_members\n")
    for i, (rep, members) in enumerate(cluster_items, start=1):
        vid = f"vOTU_{i:05d}"
        fo.write(f"{vid}\t{rep}\t{len(members)}\n")

with open(r"$TMPDIR/vOTU_rep_ids.txt", "w") as fo:
    for rep, members in cluster_items:
        fo.write(rep + "\n")
PY

run_env "$PY_ENV" python - <<PY
from pathlib import Path

keep_ids = set(x.strip() for x in open(r"$TMPDIR/vOTU_rep_ids.txt") if x.strip())
infa = Path(r"$OUTDIR/final_viral_hq.fa")
outfa = Path(r"$OUTDIR/vOTU_rep.fa")
stats = Path(r"$OUTDIR/vOTU_rep.stats.tsv")

def fasta_iter(fp):
    sid = None
    seq = []
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(seq)
                sid = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if sid is not None:
            yield sid, "".join(seq)

n = 0
total_bp = 0
with open(outfa, "w") as fo:
    for sid, seq in fasta_iter(infa):
        if sid in keep_ids:
            fo.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i+80] + "\n")
            n += 1
            total_bp += len(seq)

with open(stats, "w") as so:
    so.write("file\tnseq\ttotal_bp\n")
    so.write(f"{outfa.name}\t{n}\t{total_bp}\n")
PY

############################################################
# 9. 建丰度矩阵
############################################################
echo "[6/8] quantifying vOTU abundance with coverm ..."

ABUND_DIR="$OUTDIR/abundance_per_sample"
mkdir -p "$ABUND_DIR"

for r1 in "$READS_DIR"/*_1.fastq "$READS_DIR"/*_1.fastq.gz; do
  [[ -e "$r1" ]] || continue
  base=$(basename "$r1")
  sample=${base%_1.fastq}
  sample=${sample%_1.fastq.gz}

  r2="$READS_DIR/${sample}_2.fastq"
  [[ -s "$r2" ]] || r2="$READS_DIR/${sample}_2.fastq.gz"
  [[ -s "$r2" ]] || { echo "WARN: missing R2 for $sample"; continue; }

  run_env "$COVERM_ENV" coverm contig \
    --coupled "$r1" "$r2" \
    --reference "$OUTDIR/vOTU_rep.fa" \
    --methods tpm mean covered_fraction \
    --min-read-percent-identity 95 \
    --min-covered-fraction 0.75 \
    -t "$MAP_THREADS" \
    > "$ABUND_DIR/${sample}.tsv"
done

run_env "$PY_ENV" python - <<PY
import os, glob, re
import pandas as pd

abund_dir = r"$ABUND_DIR"
files = sorted(glob.glob(os.path.join(abund_dir, "*.tsv")))
if not files:
    raise SystemExit("No abundance files found.")

def norm(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())

def find_col(cols, target):
    nt = norm(target)
    for c in cols:
        if nt == norm(c):
            return c
    for c in cols:
        if nt in norm(c):
            return c
    return None

tpm_list = []
mean_list = []
cov_list = []

for f in files:
    sample = os.path.basename(f).rsplit(".", 1)[0]
    df = pd.read_csv(f, sep="\t")
    idcol = df.columns[0]

    tpm_col = find_col(df.columns, "tpm")
    mean_col = find_col(df.columns, "mean")
    cov_col = find_col(df.columns, "covered_fraction")

    if tpm_col is None:
        raise SystemExit(f"TPM column missing in {f}")

    tpm_list.append(df.set_index(idcol)[tpm_col].rename(sample))

    if mean_col is not None:
        mean_list.append(df.set_index(idcol)[mean_col].rename(sample))

    if cov_col is not None:
        cov_list.append(df.set_index(idcol)[cov_col].rename(sample))

tpm = pd.concat(tpm_list, axis=1).fillna(0)
tpm.index.name = "representative_id"
tpm.to_csv(r"$OUTDIR/vOTU_abundance_TPM.tsv", sep="\t")

if mean_list:
    mean = pd.concat(mean_list, axis=1).fillna(0)
    mean.index.name = "representative_id"
    mean.to_csv(r"$OUTDIR/vOTU_abundance_Mean.tsv", sep="\t")

if cov_list:
    cov = pd.concat(cov_list, axis=1).fillna(0)
    cov.index.name = "representative_id"
    cov.to_csv(r"$OUTDIR/vOTU_abundance_CoveredFraction.tsv", sep="\t")
PY

############################################################
# 10. 输出
############################################################
echo "[7/8] done"
echo "Key outputs:"
echo "  $OUTDIR/contig_tool_support.tsv"
echo "  $OUTDIR/merged_candidates.precheck.fa"
echo "  $OUTDIR/checkv_final/"
echo "  $OUTDIR/checkv_kept.tsv"
echo "  $OUTDIR/final_viral_hq.fa"
echo "  $OUTDIR/vOTU_edges.tsv"
echo "  $OUTDIR/vOTU_members.tsv"
echo "  $OUTDIR/vOTU_summary.tsv"
echo "  $OUTDIR/vOTU_rep.fa"
echo "  $OUTDIR/vOTU_abundance_TPM.tsv"
echo "  $OUTDIR/vOTU_abundance_Mean.tsv"
echo "  $OUTDIR/vOTU_abundance_CoveredFraction.tsv"


=====================================================================================================================
！！！！！！！！！！！！！组装分析流程！！！！！！！！！！！！！！！
## 1 基因预测、去冗余和定量Gene prediction, cluster & quantitfy
0107 cpu07 cssdmetahit
    wd=/home/tang/cssd
	cpu=76
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
 	# 添加分析所需的软件、脚本至环境变量，添加至~/.bashrc中自动加载，其中EasyMicrobiome文件夹放在根目录！！！
    PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/EasyMicrobiome/linux:/EasyMicrobiome/script
    echo $PATH
### metaProdigal基因预测Gene prediction
    conda activate megahit
all_samples_contigs.fasta
    mkdir -p temp/prodigal
    time prodigal -i result/megahit/all_samples_contigs.fasta \
        -d temp/prodigal/gene.fa \
        -o temp/prodigal/gene.gff \
        -p meta -f gff > temp/prodigal/gene.log 2>&1 
    # 统计基因数量,6G18s3M
    seqkit stat temp/prodigal/gene.fa 
    # 统计完整基因数量，数据量大可只用完整基因部分
    grep -c 'partial=00' temp/prodigal/gene.fa 
    # 提取完整基因(完整片段获得的基因全为完整，如成环的细菌基因组)
    grep 'partial=00' temp/prodigal/gene.fa | cut -f1 -d ' '| sed 's/>//' > temp/prodigal/full_length.id
    seqkit grep -f temp/prodigal/full_length.id temp/prodigal/gene.fa > temp/prodigal/full_length.fa
    seqkit stat temp/prodigal/full_length.fa
0107 cpu07 cssdmetahit end
### cd-hit基因聚类/去冗余cluster & redundancy
0112
    mkdir -p result/NR
    cd-hit-est -i temp/prodigal/gene.fa \
        -o result/NR/nucleotide.fa \
        -aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0
    # 统计非冗余基因数量，单次拼接结果数量下降不大，如3M-2M，多批拼接冗余度高
    grep -c '>' result/NR/nucleotide.fa
    # 翻译核酸为对应蛋白序列, --trim去除结尾的*
    seqkit translate --trim result/NR/nucleotide.fa \
        > result/NR/protein.fa 
	
#如果序列出问题执行这一步操作，对序列过滤，然后将后续protein.fa都修改为final_protein.fa！！！！！！！！！！！！！！
# 核心：仅用seqkit过滤空序列（-g）和长度<10的序列（-m 10），直接生成最终文件
seqkit seq -g -m 10 -w 0 --quiet result/NR/protein.fa > result/NR/final_protein.fa
# 1. 检查文件大小（不为0则成功）
ls -lh result/NR/final_protein.fa
# 预期输出：类似 "-rw-rw-r-- 1 aczhxzz9wi aczhxzz9wi 2.7G Dec  2 11:30 final_protein.fa"
# 2. 检查序列统计（min_len≥10，num_seqs≠0）
seqkit stats -a result/NR/final_protein.fa
# 预期输出：num_seqs 1800万+，min_len≥10，sum_len≠0

### salmon基因定量quantitfy
    mkdir -p temp/salmon
    salmon -v # 1.8.0

    # 建索引, -t序列, -i 索引，10s
    salmon index -t result/NR/nucleotide.fa \
      -p ${cpu} -i temp/salmon/index 
cpu1=38
    time tail -n+2 result/metadata.txt | cut -f1 | rush -j 2 \
      "salmon quant -i temp/salmon/index -l A -p ${cpu1} --meta \
        -1 temp/hr/{1}_1.fastq -2 temp/hr/{1}_2.fastq \
        -o temp/salmon/{1}.quant"

    # 合并
    mkdir -p result/salmon
    salmon quantmerge --quants temp/salmon/*.quant \
        -o result/salmon/gene.TPM
    salmon quantmerge --quants temp/salmon/*.quant \
        --column NumReads -o result/salmon/gene.count
    sed -i '1 s/.quant//g' result/salmon/gene.*

    # 预览结果表格
    head -n3 result/salmon/gene.*
0112 enddone---------------------------
## 2 功能基因注释Functional gene annotation

    # 输入数据：上一步预测的蛋白序列 result/NR/final_protein.fa
    # 中间结果：temp/eggnog/protein.emapper.seed_orthologs
    #           temp/eggnog/output.emapper.annotations
    #           temp/eggnog/output

    # COG定量表：result/eggnog/cogtab.count
    #            result/eggnog/cogtab.count.spf (用于STAMP)

    # KO定量表：result/eggnog/kotab.count
    #           result/eggnog/kotab.count.spf  (用于STAMP)

    # CAZy碳水化合物注释和定量：result/dbcan3/cazytab.count
    #                           result/dbcan3/cazytab.count.spf (用于STAMP)

    # 抗生素抗性：result/resfam/resfam.count
    #             result/resfam/resfam.count.spf (用于STAMP)

    # 这部分可以拓展到其它数据库

### eggNOG基因注释gene annotation(COG/KEGG/CAZy)
0116 cpu07
    # 运行并记录软件版本
    conda activate eggnog
    emapper.py --version

    # 运行emapper，18m，默认diamond 1e-3; 2M,32p,1.5h
    mkdir -p temp/eggnog
    time emapper.py --data_dir ${db}/eggnog \
      -i result/NR/final_protein.fa --cpu ${cpu} -m diamond --override \
      -o temp/eggnog/output

    # 格式化结果并显示表头
    grep -v '^##' temp/eggnog/output.emapper.annotations | sed '1 s/^#//' \
      > temp/eggnog/output
    csvtk -t headers -v temp/eggnog/output

    # 生成COG/KO/CAZy丰度汇总表
    mkdir -p result/eggnog
    # 显示帮助
    summarizeAbundance.py -h
    # 汇总，7列COG_category按字母分隔，12列KEGG_ko和19列CAZy按逗号分隔，原始值累加
    summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m temp/eggnog/output \
      -c '7,12,19' -s '*+,+,' -n raw \
      -o result/eggnog/eggnog
    sed -i 's#^ko:##' result/eggnog/eggnog.KEGG_ko.raw.txt
    sed -i '/^-/d' result/eggnog/eggnog*
    head -n3 result/eggnog/eggnog*
    # eggnog.CAZy.raw.txt  eggnog.COG_category.raw.txt  eggnog.KEGG_ko.raw.txt

    # 添加注释生成STAMP的spf格式
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      ~/EasyMicrobiome/kegg/KO_description.txt \
      result/eggnog/eggnog.KEGG_ko.raw.txt | \
      sed 's/^\t/Unannotated\t/' \
      > result/eggnog/eggnog.KEGG_ko.TPM.spf
    head -n 5 result/eggnog/eggnog.KEGG_ko.TPM.spf
    # KO to level 1/2/3
    summarizeAbundance.py \
      -i result/eggnog/eggnog.KEGG_ko.raw.txt \
      -m ~/EasyMicrobiome/kegg/KO1-4.txt \
      -c 2,3,4 -s ',+,+,' -n raw \
      -o result/eggnog/KEGG
    head -n3 result/eggnog/KEGG*
    
    # CAZy
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
       ~/EasyMicrobiome/dbcan2/CAZy_description.txt result/eggnog/eggnog.CAZy.raw.txt | \
      sed 's/^\t/Unannotated\t/' > result/eggnog/eggnog.CAZy.TPM.spf
    head -n 3 result/eggnog/eggnog.CAZy.TPM.spf
    
    # COG
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
      ~/EasyMicrobiome/eggnog/COG.anno result/eggnog/eggnog.COG_category.raw.txt > \
      result/eggnog/eggnog.COG_category.TPM.spf
    head -n 3 result/eggnog/eggnog.COG_category.TPM.spf


### CAZy碳水化合物酶/CAZyDB

    mkdir -p temp/dbcan3 result/dbcan3
    # --sensitive慢10倍，dbcan3e值为1e-102，此处以1e-3演示
    time diamond blastp \
      --db ${db}/dbcan3/CAZyDB \
      --query result/NR/final_protein.fa \
      --threads ${cpu} -e 1e-102 --outfmt 6 --max-target-seqs 1 --quiet \
      --out temp/dbcan3/gene_diamond.f6
    wc -l temp/dbcan3/gene_diamond.f6
    # 提取基因与dbcan分类对应表，按Evalue值过滤，推荐1e-102，此处演示1e-3为了有足够结果
    format_dbcan2list.pl \
      -i temp/dbcan3/gene_diamond.f6 \
      -o temp/dbcan3/gene.list 
    # 按对应表累计丰度，依赖
    summarizeAbundance.py \
      -i $wd/result/salmon/gene.TPM \
      -m $wd/temp/dbcan3/gene.list \
      -c 2 -s ',' -n raw \
      -o $wd/result/dbcan3/TPM
    # 添加注释生成STAMP的spf格式，结合metadata.txt进行差异比较
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
       ~/EasyMicrobiome/dbcan2/CAZy_description.txt result/dbcan3/TPM.CAZy.raw.txt | \
      sed 's/^\t/Unannotated\t/' \
      > result/dbcan3/TPM.CAZy.raw.spf
    head result/dbcan3/TPM.CAZy.raw.spf
    # 检查未注释数量，有则需要检查原因
    grep 'Unannotated' result/dbcan3/TPM.CAZy.raw.spf|wc -l
0112end——————————————————————————————————————————————
20260330 cpu07
    wd=/home/tang/cssd
	cpu=70
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
	fatype=all_samples_contigs.fasta #fatype=final.contigs.fa
 	# 添加分析所需的软件、脚本至环境变量，添加至~/.bashrc中自动加载，其中EasyMicrobiome文件夹放在根目录！！！
    PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/EasyMicrobiome/linux:/EasyMicrobiome/script
    echo $PATH

## 3.9 VFDB2.0建库-blast建库
# 创建ID已修复的副本
conda activate eggnog
#比对
mkdir -p temp/VFDB result/VFDB
time tblastn -query ${wd}/result/NR/final_protein.fa \
       -db ${db}/VFDB/VFDB \
       -out ${wd}/temp/VFDB/VFDB_diamond.f6 \
       -outfmt 6 \
       -evalue 1e-10 \
       -max_target_seqs 5 \
       -num_threads ${cpu}

Rscript ~/EasyMicrobiome/script/merge_diamond_abundance.R -g ${wd}/result/salmon/gene.TPM \
  -d VFDB \
  -w ${wd}/temp/VFDB \
  -o VFDB
## 3.9 VFDB2.0建库-blast建库 end---------------------------

## 3.10 SARG-diamond
#比对
mkdir -p temp/SARG result/SARG
    time diamond blastp \
      --db ${db}/SARG/SARG \
      --query ${wd}/result/NR/final_protein.fa \
      --threads ${cpu} -e 1e-10 --outfmt 6 --max-target-seqs 1 --quiet \
      --out ${wd}/temp/SARG/SARG_diamond.f6


time diamond blastp \
      --db $db/SARG/SARG \
      --query bin99.1_re.2_re.fa \
      --threads 12 -e 1e-10 --outfmt 6 --max-target-seqs 1 --quiet \
      --out ~/binningtest/bin99.1_re.2_re_diamond.f6

Rscript ~EasyMicrobiome/script/merge_diamond_abundance.R -g ${wd}/result/salmon/gene.TPM \
  -d SARG \
  -w ${wd}/temp/SARG \
  -o result/SARG
20260330 cpu07 end---------------------------
20260331 CPU01
    wd=/home/tang/cssd
	cpu=70
    soft=/home/tang/miniconda3
    db=/home/tang/db
    mkdir -p $wd && cd $wd
	fatype=all_samples_contigs.fasta #fatype=final.contigs.fa
 	# 添加分析所需的软件、脚本至环境变量，添加至~/.bashrc中自动加载，其中EasyMicrobiome文件夹放在根目录！！！
    PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/EasyMicrobiome/linux:/EasyMicrobiome/script
    echo $PATH
## 3.10 SARG-diamond
#比对
conda activate eggnog
mkdir -p temp/SARG result/SARG
    time diamond blastp \
      --db ${db}/SARG/SARG \
      --query ${wd}/result/NR/final_protein.fa \
      --threads ${cpu} -e 1e-10 --outfmt 6 --max-target-seqs 1 --quiet \
      --out ${wd}/temp/SARG/SARG_diamond.f6

Rscript ~/EasyMicrobiome/script/merge_diamond_abundance.R -g ${wd}/result/salmon/gene.TPM \
  -d SARG \
  -w ${wd}/temp/SARG \
  -o ./

# 3.11 MGE-blast  
conda activate eggnog
mkdir -p ${wd}/temp/MGE
time tblastn -query ${wd}/result/NR/final_protein.fa \
       -db ${db}/MGE/MGE \
       -out ${wd}/temp/MGE/MGE_diamond.f6 \
       -outfmt 6 \
       -evalue 1e-10 \
       -max_target_seqs 5 \
       -num_threads ${cpu}  
  
  
20260331 CPU01 end---------------------------

20260331 15:57 cctk 测试
cctk quickrun -i /home/tang/cssd/result/megahit -o /home/tang/cssd/result/cctk/cctk_quickrun_allcontig
print：
running cctk minced on the provided assemblies...
You have entered a directory name. An input file name is required: /home/tang/cssd/result/megahit/single_megahit
Total unique spacers: 6916
Total unique arrays: 766
Found 0 clusters with between 3 and 15 arrays.
Exiting.
20260402 15:57 cctk 测试 done

/home/tang/cssd/result/megahit
/home/tang/cssd/temp/bin_refinementtest/BJ/metawrap_50_10_bins
20260402 cpu01
conda activate spacepharer
#安装数据库
spacepharer downloaddb \
GenBank_phage_2018_09 \
${db}/spacepharer_phage_db \
${db}/tmp_spacepharer \
--threads 60
#再运行
mkdir -p /home/tang/cssd/result/spacepharer_out
time spacepharer easy-predict \
/home/tang/cssd/result/cctk/cctk_quickrun_allcontig/PROCESSED/CRISPR_spacers.fna \
${db}/spacepharer_phage_db \
/home/tang/cssd/result/spacepharer_out/matches.tsv \
/home/tang/cssd/result/spacepharer_out/tmp \
--threads 60 \
--fdr 0.05

20260402 cpu01 done

# 五、分箱基因组分析 =====================================
wd=/home/tang/cssd
cpu=76
cpu1=38
soft=/home/tang/miniconda3
db=/home/tang/db

mkdir -p ${wd} && cd ${wd}
mkdir -p seq temp result

fatype=all_samples_contigs.fasta   # fatype=final.contigs.fa

# 添加分析所需的软件、脚本至环境变量
PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$db/EasyMicrobiome/linux:$db/EasyMicrobiome/script
export PATH
echo $PATH



# 五、分箱基因组分析 =====================================
wd=/home/tang/cssd
cpu=76
cpu1=38
soft=/home/tang/miniconda3
db=/home/tang/db

mkdir -p ${wd} && cd ${wd}
mkdir -p seq temp result

fatype=all_samples_contigs.fasta   # fatype=final.contigs.fa

# 添加分析所需的软件、脚本至环境变量
PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$db/EasyMicrobiome/linux:$db/EasyMicrobiome/script
export PATH
echo $PATH


########################################################
## 5.1 汇总所有样本 refinement 后的 bins 到同一目录
########################################################

# 你的单样本 refinement 结果目录
bin_ref_dir=${wd}/temp/bin_refinement

# 汇总目录
mkdir -p temp/drep_in
rm -f temp/drep_in/*

# 从每个样本的 metawrap_50_10_bins 中复制 bin，并加样本名前缀避免重名
for sample_dir in ${bin_ref_dir}/*; do
    [ -d "$sample_dir" ] || continue
    sample=$(basename "$sample_dir")

    # 优先处理 .fa
    for bin in ${sample_dir}/metawrap_50_10_bins/*.fa; do
        [ -e "$bin" ] || continue
        bin_base=$(basename "$bin")
        cp "$bin" temp/drep_in/${sample}_${bin_base}
    done

    # 如果是 .fasta 也一起处理
    for bin in ${sample_dir}/metawrap_50_10_bins/*.fasta; do
        [ -e "$bin" ] || continue
        bin_base=$(basename "$bin")
        cp "$bin" temp/drep_in/${sample}_${bin_base}
    done
done

# 检查汇总结果
echo "Number of bins collected:"
ls temp/drep_in | wc -l
ls temp/drep_in | head

########################################################
## 5.2 dRep去冗余（种水平）
########################################################
20260410 cpu07
conda activate drep
cd ${wd}

mkdir -p temp/drep95

# 如果 dRep 以前跑过且想重跑，可先删除中间checkM结果
# /bin/rm -rf temp/drep95/data/checkM

# 自动判断输入文件扩展名
n_fa=$(ls temp/drep_in/*.fa 2>/dev/null | wc -l)
n_fasta=$(ls temp/drep_in/*.fasta 2>/dev/null | wc -l)

if [ "$n_fa" -gt 0 ]; then
    time dRep dereplicate temp/drep95/ \
      -g temp/drep_in/*.fa \
      -sa 0.95 -nc 0.30 -comp 50 -con 10 -p ${cpu}
elif [ "$n_fasta" -gt 0 ]; then
    time dRep dereplicate temp/drep95/ \
      -g temp/drep_in/*.fasta \
      -sa 0.95 -nc 0.30 -comp 50 -con 10 -p ${cpu}
else
    echo "Error: No genome files (.fa or .fasta) found in temp/drep_in/"
    exit 1
fi

# 查看去冗余后代表基因组数量
echo "Dereplicated genomes:"
ls temp/drep95/dereplicated_genomes/ | wc -l
ls temp/drep95/dereplicated_genomes/ | head

# 提取代表基因组ID
mkdir -p temp/drep95/data_tables
ls temp/drep95/dereplicated_genomes/ | sed 's/\.fa$//' | sed 's/\.fasta$//' \
  > temp/drep95/data_tables/dereplicated_genomes.id

# 如果 format_drep2cluster.pl 不在 PATH 中，可改成绝对路径
# 先查找该脚本位置
# find ${db} -name "format_drep2cluster.pl"

# 假设脚本位于 EasyMicrobiome/script 下，可用如下命令：
if [ -f "${db}/EasyMicrobiome/script/format_drep2cluster.pl" ]; then
    perl ${db}/EasyMicrobiome/script/format_drep2cluster.pl \
      -i temp/drep95/data_tables/Cdb.csv \
      -d temp/drep95/data_tables/dereplicated_genomes.id \
      -o temp/drep95/data_tables/Cdb.list \
      -h header num
else
    echo "Warning: format_drep2cluster.pl not found, skip Cdb.list generation"
fi

# 主要结果：
# temp/drep95/dereplicated_genomes/*.fa 或 *.fasta
# temp/drep95/data_tables/Cdb.csv
# temp/drep95/figures/*clustering*


########################################################
## 5.3 CoverM基因组定量
########################################################

conda activate coverm
mkdir -p temp/coverm
mkdir -p result/coverm

# 自动判断 dRep 输出扩展名
n_drep_fa=$(ls temp/drep95/dereplicated_genomes/*.fa 2>/dev/null | wc -l)
n_drep_fasta=$(ls temp/drep95/dereplicated_genomes/*.fasta 2>/dev/null | wc -l)

if [ "$n_drep_fa" -gt 0 ]; then
    genome_ext=fa
elif [ "$n_drep_fasta" -gt 0 ]; then
    genome_ext=fasta
else
    echo "Error: No dereplicated genomes found in temp/drep95/dereplicated_genomes/"
    exit 1
fi

# 并行计算丰度
tail -n+2 result/metadata.txt | cut -f1 | rush -j 2 \
  "coverm genome --coupled temp/hr/{}_1.fastq temp/hr/{}_2.fastq \
  -t ${cpu1} \
  --genome-fasta-directory temp/drep95/dereplicated_genomes/ \
  -x ${genome_ext} \
  -o temp/coverm/{}.txt > temp/coverm/{}.log "

# 合并结果
conda activate humann3

# 兼容 CoverM 输出表头
sed -i 's/_1.fastq Relative Abundance (%)//' temp/coverm/*.txt

humann_join_tables --input temp/coverm \
  --file_name txt \
  --output result/coverm/abundance.tsv

csvtk -t stat result/coverm/abundance.tsv

# 按组求均值（可选，没跑）
# 注意：metadata 中必须包含样本名列和 Group 列
Rscript ${db}/EasyMicrobiome/script/otu_mean.R \
  --input result/coverm/abundance.tsv \
  --metadata result/metadata.txt \
  --group Group --thre 0 \
  --scale TRUE --zoom 100 --all TRUE --type mean \
  --output result/coverm/group_mean.txt


########################################################
## 5.4 GTDB-Tk物种注释和进化树
########################################################

conda activate gtdbtk
export GTDBTK_DATA_PATH="${db}/gtdb/release226"
gtdbtk -v    # 2.3.2

mkdir -p temp/gtdb_classify

time gtdbtk classify_wf \
    --genome_dir temp/drep95/dereplicated_genomes \
    --out_dir temp/gtdb_classify \
    --extension fa \
    --skip_ani_screen \
    --prefix tax \
    --cpus ${cpu}

# 查看结果
less -S temp/gtdb_classify/tax.bac120.summary.tsv
less -S temp/gtdb_classify/tax.ar53.summary.tsv

20260413 cpu02 
########################################################
## 0. 前置准备
########################################################
wd=/home/tang/cssd
cpu=76
cpu1=38
soft=/home/tang/miniconda3
db=/home/tang/db

mkdir -p ${wd} && cd ${wd}
mkdir -p seq temp result

fatype=all_samples_contigs.fasta   # fatype=final.contigs.fa

# 添加分析所需的软件、脚本至环境变量
PATH=$soft/bin:$soft/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$db/EasyMicrobiome/linux:$db/EasyMicrobiome/script
export PATH
echo $PATH

########################################################
## 5.5 MetaCHIP 预测分箱后MAG间HGT
## 说明：
## 1) 输入基因组目录：temp/drep95/dereplicated_genomes
## 2) GTDB-Tk 分类结果目录：temp/gtdb_classify
## 3) 所有输出统一到：temp/metachip
########################################################

conda activate metachip_env
mkdir -p ${wd}/temp/metachip
cd ${wd}

# 1. 从 GTDB-Tk 结果中提取 MetaCHIP 所需 taxonomy 文件
# 仅保留前两列：user_genome 和 classification
rm -f ${wd}/temp/metachip/gtdb_taxonomy.tsv

awk 'NR==1{next} {print $1"\t"$2}' ${wd}/temp/gtdb_classify/tax.bac120.summary.tsv > ${wd}/temp/metachip/gtdb_taxonomy.tsv
awk 'NR==1{next} {print $1"\t"$2}' ${wd}/temp/gtdb_classify/tax.ar53.summary.tsv >> ${wd}/temp/metachip/gtdb_taxonomy.tsv

# 2. 简单检查 taxonomy 文件
head -n 5 ${wd}/temp/metachip/gtdb_taxonomy.tsv
wc -l ${wd}/temp/metachip/gtdb_taxonomy.tsv

# 3. 检查 taxonomy 中基因组名是否和 fasta 文件名匹配
comm -23 \
  <(cut -f1 ${wd}/temp/metachip/gtdb_taxonomy.tsv | sort) \
  <(for i in ${wd}/temp/drep95/dereplicated_genomes/*.fa; do basename "$i" .fa; done | sort)
# 如果上面没有输出，说明名字匹配，可继续运行

# 4. 运行 MetaCHIP PI
runid=mc_$(date +%Y%m%d_%H%M%S)
prefix=${runid}
parent=${wd}/temp/metachip
outdir=${parent}/${runid}

mkdir -p "${parent}"

time MetaCHIP PI \
  -p "${prefix}" \
  -r pcofg \
  -t ${cpu} \
  -o "${outdir}" \
  -i ${wd}/temp/drep95/dereplicated_genomes \
  -x fa \
  -taxon ${wd}/temp/metachip/gtdb_taxonomy.tsv \
  2>&1 | tee "${parent}/${runid}.PI.log"

#20260413done

time MetaCHIP BP \
  -p "${prefix}" \
  -r pcofg \
  -t ${cpu} \
  -o "${outdir}" \
  2>&1 | tee "${parent}/${runid}.BP.log"


# 6. 查看结果
find "${outdir}" -maxdepth 2 | sort | sed -n '1,120p'
find "${outdir}" -name "*detected_HGTs.txt"
tail -n 30 "${parent}/${runid}.PI.log"


====================================

解读consensus_MAG_ARG_summary.tsv
| 列名                       | 含义                        | 怎么理解                                                                     |
| ------------------------ | ------------------------- | ------------------------------------------------------------------------ |
| `MAG_ID`                 | MAG 的编号                   | 例如 `CD1_bin.5`，代表一个分箱基因组                                                 |
| `SARG_total`             | SARG 注释到的 ARG 基因总数        | 这个 MAG 中被 SARG 识别为 ARG 的基因数量                                             |
| `SARG_type_richness`     | SARG 注释到的 ARG 大类数量        | 这个 MAG 携带多少种抗性类型，例如 multidrug、tetracycline、aminoglycoside                |
| `SARG_subtype_richness`  | SARG 注释到的 ARG 亚型数量        | 这个 MAG 携带多少种具体 ARG subtype，例如 `tetracycline__tetA(48)`、`multidrug__mexB` |
| `DeepARG_total`          | DeepARG 注释到的 ARG 基因总数     | 这个 MAG 中被 DeepARG 识别为 ARG 的基因数量                                          |
| `DeepARG_class_richness` | DeepARG 注释到的 ARG class 数量 | 这个 MAG 携带多少种 DeepARG 抗性类别，例如 beta-lactam、multidrug、tetracycline          |


# ============================================================
# 6. MAG annotation after Prodigal
# ============================================================

wd=/home/tang/cssd
db=/home/tang/db
cpu=78

# 如果脚本中 conda activate 不生效，取消下一行注释
# source /home/tang/miniconda3/etc/profile.d/conda.sh


# ============================================================
# 6.1 Merge Prodigal faa / ffn / gff
# 输入：
#   ${wd}/temp/bin_prodigal/*.faa
#   ${wd}/temp/bin_prodigal/*.ffn
#   ${wd}/temp/bin_prodigal/*.gff
# 输出：
#   ${wd}/temp/bin_prodigal/all_MAGs.faa
#   ${wd}/temp/bin_prodigal/all_MAGs.ffn
#   ${wd}/temp/bin_prodigal/all_MAGs.gff
# ============================================================

mkdir -p ${wd}/temp/bin_prodigal

rm -f ${wd}/temp/bin_prodigal/all_MAGs.faa
rm -f ${wd}/temp/bin_prodigal/all_MAGs.ffn
rm -f ${wd}/temp/bin_prodigal/all_MAGs.gff

# merge faa
for faa in ${wd}/temp/bin_prodigal/*.faa; do
    [ "$(basename "$faa")" = "all_MAGs.faa" ] && continue
    mag=$(basename "$faa" .faa)

    awk -v m="$mag" '
    /^>/{
        sub(/^>/, ">"m"|", $0)
        print
        next
    }
    {print}
    ' "$faa" >> ${wd}/temp/bin_prodigal/all_MAGs.faa
done

# merge ffn
for ffn in ${wd}/temp/bin_prodigal/*.ffn; do
    [ "$(basename "$ffn")" = "all_MAGs.ffn" ] && continue
    mag=$(basename "$ffn" .ffn)

    awk -v m="$mag" '
    /^>/{
        sub(/^>/, ">"m"|", $0)
        print
        next
    }
    {print}
    ' "$ffn" >> ${wd}/temp/bin_prodigal/all_MAGs.ffn
done

# merge gff
echo -e "MAG_ID\tseqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes" \
> ${wd}/temp/bin_prodigal/all_MAGs.gff

for gff in ${wd}/temp/bin_prodigal/*.gff; do
    [ "$(basename "$gff")" = "all_MAGs.gff" ] && continue
    mag=$(basename "$gff" .gff)

    awk -v m="$mag" 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=9 {print m,$1,$2,$3,$4,$5,$6,$7,$8,$9}' "$gff" \
    >> ${wd}/temp/bin_prodigal/all_MAGs.gff
done

echo "===== merged Prodigal files ====="
grep -c "^>" ${wd}/temp/bin_prodigal/all_MAGs.faa
grep -c "^>" ${wd}/temp/bin_prodigal/all_MAGs.ffn
head -n 4 ${wd}/temp/bin_prodigal/all_MAGs.faa
head -n 4 ${wd}/temp/bin_prodigal/all_MAGs.ffn


# ============================================================
# 6.2 SARG annotation
# 输入：all_MAGs.faa
# 输出目录：${wd}/result/bin_SARG
# ============================================================

conda activate eggnog

mkdir -p ${wd}/result/bin_SARG

if [ ! -s ${db}/SARG/SARG_rank.txt ]; then
    echo "Error: ${db}/SARG/SARG_rank.txt not found or empty."
    echo "Please generate SARG_rank.txt first."
    exit 1
fi

time diamond blastp \
  --db ${db}/SARG/SARG \
  --query ${wd}/temp/bin_prodigal/all_MAGs.faa \
  --threads ${cpu} \
  --evalue 1e-10 \
  --max-target-seqs 1 \
  --quiet \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
  --out ${wd}/result/bin_SARG/SARG_diamond.tsv

# strict filter: identity >=70%, qcov >=70%, scov >=70%, evalue <=1e-5
awk 'BEGIN{
    FS=OFS="\t";
    print "qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore","qcov","scov"
}
{
    qcov=$4/$5
    scov=$4/$6

    if(qcov>1) qcov=1
    if(scov>1) scov=1

    if($3>=70 && qcov>=0.7 && scov>=0.7 && ($7+0)<=1e-5)
        print $1,$2,$3,$4,$5,$6,$7,$8,qcov,scov
}' ${wd}/result/bin_SARG/SARG_diamond.tsv \
> ${wd}/result/bin_SARG/SARG_diamond.filtered.tsv

# merge SARG_rank annotation
awk 'BEGIN{
    FS=OFS="\t"
}
NR==FNR{
    if(NR==1) next

    sid=$2
    type=$3
    subtype=$4
    hmm=$5
    mech_group=$6
    mech_subgroup=$7
    mech_subgroup2=$8
    rank=$12

    if(rank=="") rank="NA"

    ann[sid]=type"\t"subtype"\t"hmm"\t"mech_group"\t"mech_subgroup"\t"mech_subgroup2"\t"rank
    next
}
FNR==1{
    print "MAG_ID","Gene_ID","SARG_ID","ARG_type","ARG_subtype","HMM_category","Mechanism_group","Mechanism_subgroup","Mechanism_subgroup2","Rank","pident","length","qlen","slen","evalue","bitscore","qcov","scov"
    next
}
{
    if(!($2 in ann)) next

    split($1,a,"|")
    mag=a[1]

    print mag,$1,$2,ann[$2],$3,$4,$5,$6,$7,$8,$9,$10
}' ${db}/SARG/SARG_rank.txt \
   ${wd}/result/bin_SARG/SARG_diamond.filtered.tsv \
> ${wd}/result/bin_SARG/MAG_gene_ARG_table.tsv

# unmatched SARG IDs
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{
    if(NR==1) next
    ann[$2]=1
    next
}
FNR==1{next}
{
    if(!($2 in ann)) print $2
}' ${db}/SARG/SARG_rank.txt \
   ${wd}/result/bin_SARG/SARG_diamond.filtered.tsv \
| sort -u \
> ${wd}/result/bin_SARG/SARG_unmatched_ids.txt

# SARG subtype / type / rank / burden
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{next}
{
    mag=$1
    type=$4
    subtype=$5
    rank=$10

    if(type=="" || subtype=="") next
    if(rank=="") rank="NA"

    sub_cnt[mag,subtype]++
    type_cnt[mag,type]++
    rank_cnt[mag,rank]++

    total[mag]++
    type_set[mag,type]=1
    sub_set[mag,subtype]=1
    rank_set[mag,rank]=1
}
END{
    print "MAG_ID","ARG_subtype","Gene_count" > "'${wd}'/result/bin_SARG/MAG_ARG_subtype_count.tsv"
    for(k in sub_cnt){
        split(k,a,SUBSEP)
        print a[1],a[2],sub_cnt[k] > "'${wd}'/result/bin_SARG/MAG_ARG_subtype_count.tsv"
    }

    print "MAG_ID","ARG_type","Gene_count" > "'${wd}'/result/bin_SARG/MAG_ARG_type_count.tsv"
    for(k in type_cnt){
        split(k,a,SUBSEP)
        print a[1],a[2],type_cnt[k] > "'${wd}'/result/bin_SARG/MAG_ARG_type_count.tsv"
    }

    print "MAG_ID","ARG_rank","Gene_count" > "'${wd}'/result/bin_SARG/MAG_ARG_rank_count.tsv"
    for(k in rank_cnt){
        split(k,a,SUBSEP)
        print a[1],a[2],rank_cnt[k] > "'${wd}'/result/bin_SARG/MAG_ARG_rank_count.tsv"
    }

    print "MAG_ID","ARG_total","ARG_type_richness","ARG_subtype_richness","ARG_rank_richness" > "'${wd}'/result/bin_SARG/MAG_ARG_burden.tsv"
    for(m in total){
        tp=0
        st=0
        rk=0

        for(k in type_set){
            split(k,a,SUBSEP)
            if(a[1]==m) tp++
        }

        for(k in sub_set){
            split(k,a,SUBSEP)
            if(a[1]==m) st++
        }

        for(k in rank_set){
            split(k,a,SUBSEP)
            if(a[1]==m) rk++
        }

        print m,total[m],tp,st,rk > "'${wd}'/result/bin_SARG/MAG_ARG_burden.tsv"
    }
}' ${wd}/result/bin_SARG/MAG_gene_ARG_table.tsv

echo "===== SARG check ====="
for f in \
${wd}/result/bin_SARG/SARG_diamond.filtered.tsv \
${wd}/result/bin_SARG/MAG_gene_ARG_table.tsv \
${wd}/result/bin_SARG/MAG_ARG_type_count.tsv \
${wd}/result/bin_SARG/MAG_ARG_subtype_count.tsv \
${wd}/result/bin_SARG/MAG_ARG_rank_count.tsv \
${wd}/result/bin_SARG/MAG_ARG_burden.tsv \
${wd}/result/bin_SARG/SARG_unmatched_ids.txt
do
    echo "===== $f ====="
    head -n 5 "$f" | column -t -s $'\t'
    echo "Rows: $(wc -l < "$f")"
    echo
done


# ============================================================
# 6.3 DeepARG annotation
# 输入：all_MAGs.faa
# 输出目录：${wd}/result/bin_DeepARG
# ============================================================

conda activate deeparg_env

mkdir -p ${wd}/result/bin_DeepARG

deeparg predict \
  --model LS \
  --type prot \
  --input-file ${wd}/temp/bin_prodigal/all_MAGs.faa \
  --output-file ${wd}/result/bin_DeepARG/all_MAGs_deeparg \
  --data-path ${db}/DeepARG/deeparg

# MAG-gene-DeepARG table
awk 'BEGIN{
    FS=OFS="\t"
    print "MAG_ID","Gene_ID","ARG_name","ARG_class","Best_hit","Probability","Identity","Align_length","Bitscore","Evalue","Counts"
}
NR==1{next}
$5!="unclassified"{
    split($4,a,"|")
    mag=a[1]
    print mag,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12
}' ${wd}/result/bin_DeepARG/all_MAGs_deeparg.mapping.ARG \
> ${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.tsv

# strict DeepARG table: remove obvious regulators
awk 'BEGIN{FS=OFS="\t"}
NR==1{print; next}
{
    arg=toupper($3)
    hit=tolower($5)

    if(arg=="MEXT" || arg=="SMER" || arg=="KDPE") next
    if(hit ~ /regulator/ || hit ~ /transcriptional regulator/ || hit ~ /response regulator/) next

    print
}' ${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.tsv \
> ${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.strict.tsv

# DeepARG strict burden
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{next}
{
    total[$1]++
    class[$1 SUBSEP $4]=1
}
END{
    print "MAG_ID","DeepARG_total","DeepARG_class_richness"
    for(m in total){
        cl=0
        for(k in class){
            split(k,a,SUBSEP)
            if(a[1]==m) cl++
        }
        print m,total[m],cl
    }
}' ${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.strict.tsv \
> ${wd}/result/bin_DeepARG/MAG_DeepARG_burden.strict.tsv

echo "===== DeepARG check ====="
for f in \
${wd}/result/bin_DeepARG/all_MAGs_deeparg.mapping.ARG \
${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.tsv \
${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.strict.tsv \
${wd}/result/bin_DeepARG/MAG_DeepARG_burden.strict.tsv
do
    echo "===== $f ====="
    head -n 5 "$f" | column -t -s $'\t'
    echo "Rows: $(wc -l < "$f")"
    echo
done


# ============================================================
# 6.4 VFDB annotation
# 输入：all_MAGs.faa
# 输出目录：${wd}/result/bin_VFDB
# ============================================================

conda activate eggnog

mkdir -p ${wd}/temp/bin_VFDB ${wd}/result/bin_VFDB

time tblastn \
  -query ${wd}/temp/bin_prodigal/all_MAGs.faa \
  -db ${db}/VFDB/VFDB \
  -out ${wd}/temp/bin_VFDB/VFDB_tblastn.f6 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -evalue 1e-10 \
  -max_target_seqs 5 \
  -num_threads ${cpu}

# tblastn: slen is nucleotide length, so scov = length / (slen / 3)
awk 'BEGIN{
    FS=OFS="\t";
    print "qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore","qcov","scov"
}
{
    qcov=$4/$5
    scov=$4/($6/3)

    if(qcov>1) qcov=1
    if(scov>1) scov=1

    if($3>=60 && qcov>=0.7 && scov>=0.7 && $7<=1e-5)
        print $1,$2,$3,$4,$5,$6,$7,$8,qcov,scov
}' ${wd}/temp/bin_VFDB/VFDB_tblastn.f6 \
> ${wd}/result/bin_VFDB/VFDB_tblastn.filtered.tsv

# best VFDB hit per query gene, avoid sort -t tab issue
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{
    header=$0
    next
}
{
    q=$1
    bs=$8+0

    if(!(q in best_score) || bs > best_score[q]){
        best_score[q]=bs
        best_line[q]=$0
    }
}
END{
    print header
    for(q in best_line){
        print best_line[q]
    }
}' ${wd}/result/bin_VFDB/VFDB_tblastn.filtered.tsv \
> ${wd}/result/bin_VFDB/VFDB_tblastn.filtered.best.tsv

# annotation from the fasta used for makeblastdb
awk 'BEGIN{
    OFS="\t"
}
/^>/{
    h=substr($0,2)
    split(h,a," ")
    id=a[1]
    print id,h
}' ${db}/VFDB/eVFGC_fixed.fasta \
> ${wd}/result/bin_VFDB/VFDB_subject_annotation.tsv

# MAG-gene-VFDB table
awk 'BEGIN{
    FS=OFS="\t"
}
NR==FNR{
    vf[$1]=$2
    next
}
FNR==1{
    print "MAG_ID","Gene_ID","VFDB_subject","VFDB_description","pident","length","qlen","slen","evalue","bitscore","qcov","scov"
    next
}
{
    split($1,a,"|")
    mag=a[1]
    desc=($2 in vf ? vf[$2] : "NA")
    print mag,$1,$2,desc,$3,$4,$5,$6,$7,$8,$9,$10
}' ${wd}/result/bin_VFDB/VFDB_subject_annotation.tsv \
   ${wd}/result/bin_VFDB/VFDB_tblastn.filtered.best.tsv \
> ${wd}/result/bin_VFDB/MAG_gene_VFDB_table.tsv

# VFDB burden
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{next}
{
    mag=$1
    vf=$3
    total[mag]++
    vf_set[mag,vf]=1
}
END{
    print "MAG_ID","VFDB_total","VFDB_richness"
    for(m in total){
        n=0
        for(k in vf_set){
            split(k,a,SUBSEP)
            if(a[1]==m) n++
        }
        print m,total[m],n
    }
}' ${wd}/result/bin_VFDB/MAG_gene_VFDB_table.tsv \
> ${wd}/result/bin_VFDB/MAG_VFDB_burden.tsv

# VFDB subject count
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{next}
{
    mag=$1
    vf=$3
    desc=$4

    key=mag SUBSEP vf
    cnt[key]++

    if(!(key in description)){
        description[key]=desc
    }
}
END{
    print "MAG_ID","VFDB_subject","VFDB_description","Gene_count"
    for(k in cnt){
        split(k,a,SUBSEP)
        print a[1],a[2],description[k],cnt[k]
    }
}' ${wd}/result/bin_VFDB/MAG_gene_VFDB_table.tsv \
> ${wd}/result/bin_VFDB/MAG_VFDB_subject_count.tsv

echo "===== VFDB check ====="
for f in \
${wd}/result/bin_VFDB/VFDB_tblastn.filtered.tsv \
${wd}/result/bin_VFDB/VFDB_tblastn.filtered.best.tsv \
${wd}/result/bin_VFDB/VFDB_subject_annotation.tsv \
${wd}/result/bin_VFDB/MAG_gene_VFDB_table.tsv \
${wd}/result/bin_VFDB/MAG_VFDB_burden.tsv \
${wd}/result/bin_VFDB/MAG_VFDB_subject_count.tsv
do
    echo "===== $f ====="
    head -n 5 "$f" | column -t -s $'\t'
    echo "Rows: $(wc -l < "$f")"
    echo
done


# ============================================================
# 6.5 MGE annotation
# 输入：all_MAGs.ffn
# 输出目录：${wd}/result/bin_MGE
# ============================================================

conda activate eggnog

mkdir -p ${wd}/temp/bin_MGE ${wd}/result/bin_MGE

time blastn \
  -query ${wd}/temp/bin_prodigal/all_MAGs.ffn \
  -db ${db}/MGE/MGE \
  -out ${wd}/temp/bin_MGE/MGE_blastn.f6 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -evalue 1e-10 \
  -max_target_seqs 5 \
  -num_threads ${cpu}

# filter: identity >=80%, qcov >=70%, scov >=70%, evalue <=1e-5
awk 'BEGIN{
    FS=OFS="\t";
    print "qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore","qcov","scov"
}
{
    qcov=$4/$5
    scov=$4/$6

    if(qcov>1) qcov=1
    if(scov>1) scov=1

    if($3>=80 && qcov>=0.7 && scov>=0.7 && $7<=1e-5)
        print $1,$2,$3,$4,$5,$6,$7,$8,qcov,scov
}' ${wd}/temp/bin_MGE/MGE_blastn.f6 \
> ${wd}/result/bin_MGE/MGE_blastn.filtered.tsv

# best MGE hit per query gene, avoid sort -t tab issue
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{
    header=$0
    next
}
{
    q=$1
    bs=$8+0

    if(!(q in best_score) || bs > best_score[q]){
        best_score[q]=bs
        best_line[q]=$0
    }
}
END{
    print header
    for(q in best_line){
        print best_line[q]
    }
}' ${wd}/result/bin_MGE/MGE_blastn.filtered.tsv \
> ${wd}/result/bin_MGE/MGE_blastn.filtered.best.tsv

# annotation from MGE fasta
awk 'BEGIN{
    OFS="\t"
}
/^>/{
    h=substr($0,2)
    split(h,a," ")
    id=a[1]
    print id,h
}' ${db}/MGE/MGEs_FINAL_99perc_trim.fasta \
> ${wd}/result/bin_MGE/MGE_subject_annotation.tsv

# MAG-gene-MGE table
awk 'BEGIN{
    FS=OFS="\t"
}
NR==FNR{
    ann[$1]=$2
    next
}
FNR==1{
    print "MAG_ID","Gene_ID","MGE_subject","MGE_description","pident","length","qlen","slen","evalue","bitscore","qcov","scov"
    next
}
{
    split($1,a,"|")
    mag=a[1]
    desc=($2 in ann ? ann[$2] : "NA")
    print mag,$1,$2,desc,$3,$4,$5,$6,$7,$8,$9,$10
}' ${wd}/result/bin_MGE/MGE_subject_annotation.tsv \
   ${wd}/result/bin_MGE/MGE_blastn.filtered.best.tsv \
> ${wd}/result/bin_MGE/MAG_gene_MGE_table.tsv

# MGE burden
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{next}
{
    mag=$1
    mge=$3
    total[mag]++
    mge_set[mag,mge]=1
}
END{
    print "MAG_ID","MGE_total","MGE_richness"
    for(m in total){
        n=0
        for(k in mge_set){
            split(k,a,SUBSEP)
            if(a[1]==m) n++
        }
        print m,total[m],n
    }
}' ${wd}/result/bin_MGE/MAG_gene_MGE_table.tsv \
> ${wd}/result/bin_MGE/MAG_MGE_burden.tsv

# MGE subject count
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{next}
{
    mag=$1
    mge=$3
    desc=$4

    key=mag SUBSEP mge
    cnt[key]++

    if(!(key in description)){
        description[key]=desc
    }
}
END{
    print "MAG_ID","MGE_subject","MGE_description","Gene_count"
    for(k in cnt){
        split(k,a,SUBSEP)
        print a[1],a[2],description[k],cnt[k]
    }
}' ${wd}/result/bin_MGE/MAG_gene_MGE_table.tsv \
> ${wd}/result/bin_MGE/MAG_MGE_subject_count.tsv

echo "===== MGE check ====="
for f in \
${wd}/result/bin_MGE/MGE_blastn.filtered.tsv \
${wd}/result/bin_MGE/MGE_blastn.filtered.best.tsv \
${wd}/result/bin_MGE/MGE_subject_annotation.tsv \
${wd}/result/bin_MGE/MAG_gene_MGE_table.tsv \
${wd}/result/bin_MGE/MAG_MGE_burden.tsv \
${wd}/result/bin_MGE/MAG_MGE_subject_count.tsv
do
    echo "===== $f ====="
    head -n 5 "$f" | column -t -s $'\t'
    echo "Rows: $(wc -l < "$f")"
    echo
done


# ============================================================
# 6.6 MAG-level intersections
# 输出目录：${wd}/result/bin_intersect
# 先做单独注释，最后统一做 MAG 交集
# ============================================================

mkdir -p ${wd}/result/bin_intersect

# ---------- 6.6.1 SARG ∩ DeepARG strict = high-confidence ARG-host MAGs ----------

mkdir -p ${wd}/result/bin_intersect/ARG_consensus

cut -f1 ${wd}/result/bin_SARG/MAG_ARG_burden.tsv | tail -n +2 | sort -u \
> ${wd}/result/bin_intersect/ARG_consensus/SARG_MAG_ids.txt

cut -f1 ${wd}/result/bin_DeepARG/MAG_DeepARG_burden.strict.tsv | tail -n +2 | sort -u \
> ${wd}/result/bin_intersect/ARG_consensus/DeepARG_MAG_ids.txt

comm -12 \
${wd}/result/bin_intersect/ARG_consensus/SARG_MAG_ids.txt \
${wd}/result/bin_intersect/ARG_consensus/DeepARG_MAG_ids.txt \
> ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{keep[$1]=1; next}
FNR==1 || ($1 in keep)
' ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt \
  ${wd}/result/bin_SARG/MAG_ARG_burden.tsv \
> ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_burden_from_SARG.tsv

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{keep[$1]=1; next}
FNR==1 || ($1 in keep)
' ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt \
  ${wd}/result/bin_DeepARG/MAG_DeepARG_burden.strict.tsv \
> ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_burden_from_DeepARG.tsv

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{
    if(FNR>1) s[$1]=$2"\t"$3"\t"$4"\t"$5
    next
}
FNR==1{
    print "MAG_ID","SARG_total","SARG_type_richness","SARG_subtype_richness","SARG_rank_richness","DeepARG_total","DeepARG_class_richness"
    next
}
{
    if($1 in s){
        print $1,s[$1],$2,$3
    }
}' ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_burden_from_SARG.tsv \
   ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_burden_from_DeepARG.tsv \
> ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_summary.tsv

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{keep[$1]=1; next}
FNR==1 || ($1 in keep)
' ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt \
  ${wd}/result/bin_SARG/MAG_gene_ARG_table.tsv \
> ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_gene_from_SARG.tsv

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{keep[$1]=1; next}
FNR==1 || ($1 in keep)
' ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt \
  ${wd}/result/bin_DeepARG/MAG_gene_DeepARG_table.strict.tsv \
> ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_gene_from_DeepARG.tsv


# ---------- 6.6.2 Prepare VFDB and MGE MAG ID lists ----------

cut -f1 ${wd}/result/bin_VFDB/MAG_VFDB_burden.tsv | tail -n +2 | sort -u \
> ${wd}/result/bin_intersect/VFDB_MAG_ids.txt

cut -f1 ${wd}/result/bin_MGE/MAG_MGE_burden.tsv | tail -n +2 | sort -u \
> ${wd}/result/bin_intersect/MGE_MAG_ids.txt


# ---------- 6.6.3 ARG-MGE MAG intersection ----------

mkdir -p ${wd}/result/bin_intersect/ARG_MGE

comm -12 \
${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt \
${wd}/result/bin_intersect/MGE_MAG_ids.txt \
> ${wd}/result/bin_intersect/ARG_MGE/ARG_MGE_MAG_ids.txt

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{
    if(FNR>1) mge[$1]=$2"\t"$3
    next
}
FNR==1{
    print $0,"MGE_total","MGE_richness"
    next
}
{
    if($1 in mge) print $0,mge[$1]
}' ${wd}/result/bin_MGE/MAG_MGE_burden.tsv \
   ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_summary.tsv \
> ${wd}/result/bin_intersect/ARG_MGE/ARG_MGE_MAG_summary.tsv


# ---------- 6.6.4 ARG-VFDB MAG intersection ----------

mkdir -p ${wd}/result/bin_intersect/ARG_VFDB

comm -12 \
${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt \
${wd}/result/bin_intersect/VFDB_MAG_ids.txt \
> ${wd}/result/bin_intersect/ARG_VFDB/ARG_VFDB_MAG_ids.txt

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{
    if(FNR>1) vf[$1]=$2"\t"$3
    next
}
FNR==1{
    print $0,"VFDB_total","VFDB_richness"
    next
}
{
    if($1 in vf) print $0,vf[$1]
}' ${wd}/result/bin_VFDB/MAG_VFDB_burden.tsv \
   ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_summary.tsv \
> ${wd}/result/bin_intersect/ARG_VFDB/ARG_VFDB_MAG_summary.tsv


# ---------- 6.6.5 VFDB-MGE MAG intersection ----------

mkdir -p ${wd}/result/bin_intersect/VFDB_MGE

comm -12 \
${wd}/result/bin_intersect/VFDB_MAG_ids.txt \
${wd}/result/bin_intersect/MGE_MAG_ids.txt \
> ${wd}/result/bin_intersect/VFDB_MGE/VFDB_MGE_MAG_ids.txt

awk 'BEGIN{FS=OFS="\t"}
ARGIND==1{
    if(FNR>1) vf[$1]=$2"\t"$3
    next
}
ARGIND==2{
    if(FNR>1) mge[$1]=$2"\t"$3
    next
}
ARGIND==3{
    if(FNR==1){
        print "MAG_ID","VFDB_total","VFDB_richness","MGE_total","MGE_richness"
        next
    }
    if(($1 in vf) && ($1 in mge)){
        print $1,vf[$1],mge[$1]
    }
}' ${wd}/result/bin_VFDB/MAG_VFDB_burden.tsv \
   ${wd}/result/bin_MGE/MAG_MGE_burden.tsv \
   ${wd}/result/bin_intersect/VFDB_MGE/VFDB_MGE_MAG_ids.txt \
> ${wd}/result/bin_intersect/VFDB_MGE/VFDB_MGE_MAG_summary.tsv


# ---------- 6.6.6 ARG-VFDB-MGE MAG intersection ----------

mkdir -p ${wd}/result/bin_intersect/ARG_VFDB_MGE

comm -12 \
${wd}/result/bin_intersect/ARG_MGE/ARG_MGE_MAG_ids.txt \
${wd}/result/bin_intersect/VFDB_MAG_ids.txt \
> ${wd}/result/bin_intersect/ARG_VFDB_MGE/ARG_VFDB_MGE_MAG_ids.txt

awk 'BEGIN{FS=OFS="\t"}
ARGIND==1{
    if(FNR>1) vf[$1]=$2"\t"$3
    next
}
ARGIND==2{
    if(FNR>1) mge[$1]=$2"\t"$3
    next
}
ARGIND==3{
    if(FNR==1){
        print $0,"VFDB_total","VFDB_richness","MGE_total","MGE_richness"
        next
    }

    if(($1 in vf) && ($1 in mge)){
        print $0,vf[$1],mge[$1]
    }
}' ${wd}/result/bin_VFDB/MAG_VFDB_burden.tsv \
   ${wd}/result/bin_MGE/MAG_MGE_burden.tsv \
   ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_summary.tsv \
> ${wd}/result/bin_intersect/ARG_VFDB_MGE/ARG_VFDB_MGE_MAG_summary.tsv


# ============================================================
# 6.7 Final check
# ============================================================

echo "===== final MAG numbers ====="
echo -n "SARG MAGs: "
wc -l ${wd}/result/bin_intersect/ARG_consensus/SARG_MAG_ids.txt

echo -n "DeepARG strict MAGs: "
wc -l ${wd}/result/bin_intersect/ARG_consensus/DeepARG_MAG_ids.txt

echo -n "ARG consensus MAGs: "
wc -l ${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_ids.txt

echo -n "VFDB MAGs: "
wc -l ${wd}/result/bin_intersect/VFDB_MAG_ids.txt

echo -n "MGE MAGs: "
wc -l ${wd}/result/bin_intersect/MGE_MAG_ids.txt

echo -n "ARG-MGE MAGs: "
wc -l ${wd}/result/bin_intersect/ARG_MGE/ARG_MGE_MAG_ids.txt

echo -n "ARG-VFDB MAGs: "
wc -l ${wd}/result/bin_intersect/ARG_VFDB/ARG_VFDB_MAG_ids.txt

echo -n "VFDB-MGE MAGs: "
wc -l ${wd}/result/bin_intersect/VFDB_MGE/VFDB_MGE_MAG_ids.txt

echo -n "ARG-VFDB-MGE MAGs: "
wc -l ${wd}/result/bin_intersect/ARG_VFDB_MGE/ARG_VFDB_MGE_MAG_ids.txt

echo "===== final summary preview ====="
for f in \
${wd}/result/bin_intersect/ARG_consensus/ARG_consensus_MAG_summary.tsv \
${wd}/result/bin_intersect/ARG_MGE/ARG_MGE_MAG_summary.tsv \
${wd}/result/bin_intersect/ARG_VFDB/ARG_VFDB_MAG_summary.tsv \
${wd}/result/bin_intersect/VFDB_MGE/VFDB_MGE_MAG_summary.tsv \
${wd}/result/bin_intersect/ARG_VFDB_MGE/ARG_VFDB_MGE_MAG_summary.tsv
do
    echo "===== $f ====="
    head -n 5 "$f" | column -t -s $'\t'
    echo "Rows: $(wc -l < "$f")"
    echo
done

#MGA FUNCTION
20260422 报错
wd=/home/tang/cssd
cpu=76
cpu1=38
soft=/home/tang/miniconda3
db=/home/tang/db

mkdir -p ${wd} && cd ${wd}
mkdir -p seq temp result

fatype=all_samples_contigs.fasta   # fatype=final.contigs.fa

source ${soft}/etc/profile.d/conda.sh


# ============================================================
# 7. MAG functional annotation: basic tables + eggNOG
# ============================================================

mkdir -p \
  ${wd}/result/bin_MAG_function \
  ${wd}/temp/bin_eggnog \
  ${wd}/result/bin_eggnog

# 7.1 gene2MAG table
grep "^>" ${wd}/temp/bin_prodigal/all_MAGs.faa \
  | sed 's/^>//; s/ #.*//' \
  | awk 'BEGIN{FS=OFS="\t"; print "gene_id","MAG_ID"} {split($1,a,"|"); print $1,a[1]}' \
  > ${wd}/result/bin_MAG_function/gene2MAG.tsv

# 7.2 gene count per MAG
grep "^>" ${wd}/temp/bin_prodigal/all_MAGs.faa \
  | sed 's/^>//; s/ #.*//' \
  | awk 'BEGIN{FS=OFS="\t"} {split($1,a,"|"); print a[1]}' \
  | sort | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "MAG_ID","gene_count"} {print $2,$1}' \
  > ${wd}/result/bin_MAG_function/MAG_gene_count.tsv

# 7.3 GTDB-Tk taxonomy table
echo -e "MAG_ID\tGTDB_taxonomy" > ${wd}/result/bin_MAG_function/MAG_gtdb_taxonomy.tsv

for f in ${wd}/result/gtdb_classify/tax.*.summary.tsv; do
    awk 'NR>1{print $1"\t"$2}' "$f"
done | sort -u >> ${wd}/result/bin_MAG_function/MAG_gtdb_taxonomy.tsv

# 7.4 CoverM abundance table
awk 'BEGIN{FS=OFS="\t"}
NR==1{$1="MAG_ID"; print; next}
{
    g=$1;
    sub(/^.*\//,"",g);
    sub(/\.fa$/,"",g);
    sub(/\.fasta$/,"",g);
    $1=g;
    print;
}' ${wd}/result/coverm/abundance.tsv \
> ${wd}/result/bin_MAG_function/MAG_abundance.tsv

# 7.5 MAG basic information table
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{if(NR>1) gene[$1]=$2; next}
FNR==1{print "MAG_ID","GTDB_taxonomy","gene_count"; next}
{print $1,$2,gene[$1]}
' ${wd}/result/bin_MAG_function/MAG_gene_count.tsv \
  ${wd}/result/bin_MAG_function/MAG_gtdb_taxonomy.tsv \
> ${wd}/result/bin_MAG_function/MAG_basic_info.tsv

# 7.6 eggNOG-mapper annotation
conda activate eggnog

emapper.py \
  -i ${wd}/temp/bin_prodigal/all_MAGs.faa \
  --itype proteins \
  -o all_MAGs \
  --output_dir ${wd}/temp/bin_eggnog \
  --cpu ${cpu} \
  --override \
  --data_dir ${db}/eggnog

# 7.7 Add MAG_ID to eggNOG results
{
    echo -ne "MAG_ID\t"
    grep "^#query" ${wd}/temp/bin_eggnog/all_MAGs.emapper.annotations | tail -n 1 | sed 's/^#//'
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{if(NR>1) mag[$1]=$2; next}
    !/^#/ && NF>1{print mag[$1],$0}
    ' ${wd}/result/bin_MAG_function/gene2MAG.tsv \
      ${wd}/temp/bin_eggnog/all_MAGs.emapper.annotations
} > ${wd}/result/bin_eggnog/all_MAGs.emapper.annotations.withMAG.tsv

head ${wd}/result/bin_MAG_function/gene2MAG.tsv
head ${wd}/result/bin_MAG_function/MAG_gene_count.tsv
head ${wd}/result/bin_MAG_function/MAG_gtdb_taxonomy.tsv
head ${wd}/result/bin_MAG_function/MAG_abundance.tsv
head ${wd}/result/bin_MAG_function/MAG_basic_info.tsv
head ${wd}/result/bin_eggnog/all_MAGs.emapper.annotations.withMAG.tsv


# ============================================================
# 9. KofamScan KEGG KO annotation
# ============================================================

conda activate kofamscan
mkdir -p ${wd}/temp/bin_kofam ${wd}/result/bin_kofam

exec_annotation \
  -p ${db}/kofam/profiles \
  -k ${db}/kofam/ko_list \
  --cpu ${cpu} \
  -f mapper \
  -o ${wd}/temp/bin_kofam/all_MAGs.kofam.mapper.tsv \
  ${wd}/temp/bin_prodigal/all_MAGs.faa

cp ${wd}/temp/bin_kofam/all_MAGs.kofam.mapper.tsv \
   ${wd}/result/bin_kofam/

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{if(NR>1) mag[$1]=$2; next}
{print mag[$1],$1,$2}
' ${wd}/result/bin_MAG_function/gene2MAG.tsv \
  ${wd}/result/bin_kofam/all_MAGs.kofam.mapper.tsv \
> ${wd}/result/bin_kofam/all_MAGs.kofam.mapper.withMAG.tsv

sed -i '1iMAG_ID\tgene_id\tKO' \
  ${wd}/result/bin_kofam/all_MAGs.kofam.mapper.withMAG.tsv

head ${wd}/result/bin_kofam/all_MAGs.kofam.mapper.withMAG.tsv



# ============================================================
# 15. antiSMASH secondary metabolite BGC annotation
# ============================================================

conda activate antismash
antismash --version

mkdir -p ${wd}/temp/bin_antismash ${wd}/result/bin_antismash
 --databases ${db}/antismash \这里可能会报错，尝试把默认数据库复制到这个文件夹！！！！！！！！！！！！！！！！！！！！
ls ${wd}/temp/drep95/dereplicated_genomes/*.fa | ${db}/EasyMicrobiome/linux/rush -j 6 \
  "b=\$(basename {} .fa); \
   /home/tang/miniconda3/envs/antismash/bin/antismash {} \
   --databases ${db}/antismash \ 
   --output-dir ${wd}/temp/bin_antismash/\${b} \
   --genefinding-tool prodigal \
   --cpus 12 \
   --cb-general \
   --cb-knownclusters"

# BGC count per MAG
echo -e "MAG_ID\tBGC_count" > ${wd}/result/bin_antismash/MAG_BGC_count.tsv

for d in ${wd}/temp/bin_antismash/*; do
    [ -d "$d" ] || continue
    b=$(basename "$d")
    n=$(find "$d" -name "*.region*.gbk" | wc -l)
    echo -e "${b}\t${n}" >> ${wd}/result/bin_antismash/MAG_BGC_count.tsv
done

# BGC product type per MAG
echo -e "MAG_ID\tregion_file\tproduct" > ${wd}/result/bin_antismash/MAG_BGC_products.tsv

for gbk in ${wd}/temp/bin_antismash/*/*.region*.gbk; do
    mag=$(basename "$(dirname "$gbk")")
    product=$(grep -m 1 '/product=' "$gbk" | sed 's/.*\/product=//; s/"//g')
    echo -e "${mag}\t$(basename "$gbk")\t${product}" \
      >> ${wd}/result/bin_antismash/MAG_BGC_products.tsv
done

head ${wd}/result/bin_antismash/MAG_BGC_count.tsv
head ${wd}/result/bin_antismash/MAG_BGC_products.tsv
=====================================================================================================================
对198样本处理
1230
tmux new -s srr198  cpu04 220GB
cd /home/tang/lixuechen251230/198

for i in SRR* ERR*; do
    fasterq-dump -p -e 80 --split-3 -O ./ "$i"
done
#0101enddone___________________________________

# 1. 移动双端fastq文件到pair目录
find /home/tang/lixuechen251230/198 -maxdepth 1 -type f -name "*_[12].fastq" -exec mv {} /home/tang/lixuechen251230/pair/ \;

# 2. 移动单端fastq文件到singer目录
find /home/tang/lixuechen251230/198 -maxdepth 1 -type f -name "*.fastq" ! -name "*_[12].fastq" -exec mv {} /home/tang/lixuechen251230/singer/ \;

# 2. 移动单端fastq文件到singer目录
find /home/tang/lixuechen251230/singer -maxdepth 1 -type f -name "*.fastq" ! -name "*_[12].fastq" -exec mv {} /home/tang/lixuechen251230/pair/ \;

#针对单双端测序结果
0104
tmux new -s 198argoap  cpu04
wd=/home/tang/lixuechen251230
cpu=80
cd $wd
conda activate argoap
time args_oap stage_one -i $wd/pair -o $wd/argoap/output -f fastq -t $cpu  # 第一阶段，生成元数据等
time args_oap stage_two -i $wd/argoap/output -t $cpu  # 第二阶段，生成归一化结果
0104  end——————————————————————————————————

=====================================================================================================================
icamp 测试
构建系统发育树
/home/tang/cssd/temp/gtdb_classify/align/test
tax.bac120.user_msa.fasta.gz  解压
gunzip -c tax.bac120.user_msa.fasta.gz > mag_bac120_user_msa.faa
conda activate fasttree
FastTree -lg mag_bac120_user_msa.faa > mag_tree.nwk
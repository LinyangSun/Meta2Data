# JSON 参数指南

AmpliconPIP 和 AmpliconTAXA 使用 `--parameter settings.json` 加载高级参数。下表是 **JSON 字段，不是命令行选项**；例如 `primer.window` 应写成 `{"primer": {"window": 20}}`。

只需填写要修改的字段，其余使用默认值。相同设置的优先级为：**显式命令行参数 > JSON 配置 > 默认值**。配置按所选模块、平台和处理阶段生效；例如 PIP 不使用 `taxa.*`，`dada2` 不使用 `vsearch.*`。未知字段、错误类型和超出范围的值会报错。

表中范围均包含端点；“整数”字段不能填写小数，布尔值使用 `true` / `false`。**比例用 `0–1`，百分比用 `0–100`**：`0.97` 表示 97%，而接头保护的相似度应填写 `98`，不是 `0.98`。

**接头保护：`adapter_guard`**

检查 fastp 推断的候选接头是否可能是真实的 16S 序列。只有候选片段足够长、相似度和覆盖率达标，且比对 E-value ≤ `1e-5` 时，才触发更保守的重新处理。

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `adapter_guard.enabled` | `true` | 开启接头误剪检查；设为 `false` 仍执行常规去接头 | 布尔值 |
| `adapter_guard.min_length` | `50` | 可进行保护判定的候选接头最短长度，单位 bp；更短的候选不触发保护 | 整数 ≥ 1 |
| `adapter_guard.min_identity` | `98.0` | 候选接头与 16S 参考的最低比对相似度，单位 % | 数值 `0–100` |
| `adapter_guard.min_coverage` | `95.0` | 单条比对片段至少覆盖候选接头长度的百分比，单位 % | 数值 `0–100` |

`--adapter-guard` / `--no-adapter-guard` 可覆盖 JSON 中的开关。未命中参考不代表该片段一定是接头。详细流程见 [接头保护说明](adapter_guard.md)。

**自动引物检测：`primer`**

程序按样本名排序，使用第一个样本及其所有 lane/chunk 中通过检测筛选的 reads 推断引物，再将该结果应用于整个数据集。以下 `min_*` 字段只筛选**参与引物推断的 reads**，不会据此删除整个数据集中的 reads；不通过这些筛选的 reads 仍可能进入后续剪切和质控。

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `primer.window` | `20` | 检查 read 前端的碱基数；用于构建共识、匹配引物库和计算 fold，不等于固定剪切长度 | 整数 `1–100` |
| `primer.fold_threshold` | `16` | 未匹配引物库时，`fold < 此值` 判为未知引物；否则判为无引物 | 整数 ≥ 1 |
| `primer.support_frequency` | `0.10` | 每个位点上，一个碱基被视为有支持所需的最低频率；`0.10` 为 10% | 数值 `0.001–0.25` |
| `primer.database_identity` | `0.85` | 共识与引物库匹配时的最低一致比例；只计算共识中非 `N` 的有效位点，兼容 IUPAC 简并碱基 | 数值 `0–1` |
| `primer.informative_fraction` | `0.50` | 数据库匹配中，共识的有效位点数至少占对应库引物长度的比例 | 数值 `0.01–1` |
| `primer.unknown_trim_length` | `20` | 判为未知引物且继续处理时，从相应 read 前端剪除的长度，单位 bp | 整数 `1–100` |
| `primer.skip_unknown` | `false` | 判为未知引物时跳过整个数据集；可由 `--skip-unknown-primers` 开启 | 布尔值 |
| `primer.min_length` | `50` | 参与引物推断的 read 最短长度，单位 bp | 整数 ≥ 1 |
| `primer.min_average_quality` | `20.0` | 参与推断的 read 最低平均 Phred 质量；检测到单一占位质量分数时不执行此项筛选 | 数值 `0–93` |
| `primer.min_complexity` | `0.3` | 参与推断的 read 最低序列复杂度，计算为不同二碱基组合数除以 16 | 数值 `0–1` |
| `primer.min_entropy` | `1.0` | 参与推断的 read 最低 A/C/G/T 组成 Shannon 熵，单位 bit；用于排除组成过于单一的序列 | 数值 `0–2` |

每个位点支持 1、2、3、4 种碱基，分别给 fold 贡献乘数 1、2、3、4；fold 是检测窗口内这些乘数的乘积。已匹配引物库时，优先按匹配终点剪切，不再使用未知引物的 fold 判定。

这些字段只影响自动检测。为项目或本地 CSV 行手动指定引物后，程序使用指定引物，不调用自动检测。**目前没有整体跳过引物检测和剪切的用户开关**；`unknown_trim_length` 不能设为 `0`，`skip_unknown=false` 也不表示关闭引物处理。PacBio `dada2` 需要可用于定向的已知或显式引物，不能通过未知引物的固定剪切替代。

手动指定引物（项目的 `--public-primer-fwd` / `--public-primer-rev`，或本地 CSV 已映射的引物列）时，实际作用如下。反向引物的含义取决于平台和读段布局，不能一概理解为剪除每条 read 的 3′ 端。

| 数据类型 | 正向引物 | 反向引物 |
|---|---|---|
| Illumina 双端 | 用于 R1 的 5′ 端剪切（cutadapt `-g`） | 用于 R2 的 5′ 端剪切（cutadapt `-G`） |
| 普通单端：Illumina、454、Ion Torrent、ONT | 仅提供正向引物时剪切 5′ 端；同时提供两条引物时使用 linked adapter | 以其反向互补序列匹配 3′ 端；未匹配时仍可剪除正向引物 |
| PacBio `vsearch` | 同单端规则，并双方向搜索、统一 read 方向（`--revcomp`），保留原 read ID | 同单端 3′ 端规则；未提供时不自动推断或补充 |
| PacBio `dada2` | 前置步骤只记录引物，随后传给 `denoise-ccs --p-front` | 提供时传给 `denoise-ccs --p-adapter` |

单端 linked adapter 要求匹配正向引物，反向引物匹配可选。这决定是否剪切；未匹配正向引物的 reads 仍保留，不因此丢弃。

**选择 `vsearch` 时：`vsearch`**

Illumina、Ion Torrent、PacBio 使用共享的预聚类、去噪及聚类分支；454 使用下述独立顺序，仍读取共享的相似度和去噪丰度参数；ONT 使用独立聚类参数。各平台的前处理参数仅在表中注明的条件下生效。

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `vsearch.maxee` | `1.0` | Illumina 正常质量分支：每条 read 允许的最大期望错误数；双端拼接成功时作用于拼接后的 read | 数值 ≥ 0 |
| `vsearch.merge_min_fraction` | `0.5` | Illumina 正常质量双端数据：每个样本接受拼接结果的最低成功比例；低于此值则改用该样本 R1 | 数值 `0–1` |
| `vsearch.precluster_identity` | `0.99` | 共享及 454 分支进入去噪前的预聚类相似度比例 | 数值 `0.01–1` |
| `vsearch.cluster_identity` | `0.97` | 共享及 454 分支最终聚类和 reads 回贴的相似度比例 | 数值 `0.01–1` |
| `vsearch.minsize` | `2` | 共享分支预聚类前及去噪时的最低丰度；454 仅在 99% 预聚类后的去噪步骤使用 | 整数 ≥ 1 |
| `vsearch.min_frequency` | `2` | 最终特征在所有样本中合计的最低丰度，包括 ONT；454 不使用此过滤 | 整数 ≥ 1 |
| `vsearch.degraded_trim_left` | `0` | Illumina 退化、分箱质量分支：额外剪除 read 前端的长度，单位 bp | 整数 ≥ 0 |
| `vsearch.min_length` | `50` | Illumina 退化、分箱质量分支：前端剪切后保留 read 的最短长度，单位 bp | 整数 ≥ 1 |
| `vsearch.max_n` | `1` | Illumina 退化、分箱质量分支：每条保留 read 允许的最多 `N` 数量 | 整数 ≥ 0 |
| `vsearch.ion_maxee` | `null` | Ion Torrent 最大期望错误数；`null` 按数据集读长中位数和 Q19 自动计算，数值则固定覆盖 | `null` 或数值 ≥ 0 |
| `vsearch.ion_trim_left` | `0` | Ion Torrent 分支：引物处理后额外剪除前端的长度，单位 bp | 整数 ≥ 0 |
| `vsearch.pacbio_maxee_rate` | `0.01` | PacBio 分支：每碱基最大期望错误率；`0.01` 为 1%，按 read 长度计算 | 数值 `0–1` |
| `vsearch.pacbio_min_length` | `1000` | PacBio 分支的 read 长度下限，单位 bp | 整数 ≥ 1 |
| `vsearch.pacbio_max_length` | `2000` | PacBio 分支的 read 长度上限，单位 bp | 整数 ≥ 1 |

`maxee` 是根据质量分数计算的期望错误数，不是实际错配个数。PacBio 最短长度不得大于最长长度。454 不使用 Illumina 的 `maxee` 或退化质量参数。

**454 相对长度过滤：`ls454`**

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `ls454.length_fraction` | `0.5` | 每个样本的最短长度为 `ceil(引物处理及剪尾后读长中位数 × 比例)` | 数值 `0.01–1` |
| `ls454.max_n` | `1` | 每条 read 允许的最多 N 数量，合计大写 `N` 和小写 `n` | 整数 ≥ 0 |

例如 `{"ls454": {"length_fraction": 0.5, "max_n": 1}}`。中位数按每个 Run／本地样本的全部 reads 计算，包含随后可能因 N 数量超限被过滤的 reads；偶数数量取中间两个长度的平均值。先删除低于长度阈值的 reads，再删除 `N`/`n` 合计数量超过 `max_n` 的 reads，两项损失互斥。阈值取整后边界长度保留；不固定裁切保留序列、不使用 Phred 或 EE 过滤。样本内混合扩增区域或多数 reads 已经截短时，相对长度不能可靠判断真实目标范围，应结合实验信息检查结果。

过滤位于完全去重之前，逐样本输入数、中位数、阈值、两类损失、保留数及其中含 N 的数量写入 `ls454_quality-vsearch.json`。全部样本过滤后无 reads 时明确报错。默认允许含 1 个 N 的合格 reads 参与后续聚类、去噪和回贴；`max_n=0` 恢复严格去除含 N reads。

454 先对全部合格序列去重，再按 `vsearch.precluster_identity` 预聚类、累加丰度，然后执行 `vsearch.minsize` 对应的去噪、去嵌合体和最终聚类。不会在预聚类前丢弃丰度为 1 的候选，也不会因单样本独特序列比例高而跳过。回贴使用通过长度/N 过滤的 reads，并根据聚类成员关系排除明确识别的嵌合体成员；这不等于所有未识别的嵌合体都已去除。对应统计写入 `ls454_members-vsearch.json`。最终保留 count 为 1 的特征，不应用 `vsearch.min_frequency`。

**选择 `dada2` 时：`dada2`**

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `dada2.ion_maxee` | `null` | Ion Torrent：自动按读长中位数和 Q19 计算并传给 `denoise-pyro --p-max-ee`；数值则固定覆盖 | `null` 或数值 ≥ 0 |
| `dada2.ion_trim_left` | `0` | Ion Torrent：传给 `denoise-pyro` 的额外前端剪切长度，单位 bp | 整数 ≥ 0 |
| `dada2.pacbio_min_length` | `1000` | PacBio：传给 `denoise-ccs` 的最短保留长度，单位 bp | 整数 ≥ 1 |
| `dada2.pacbio_max_length` | `1600` | PacBio：传给 `denoise-ccs` 的最长保留长度，单位 bp | 整数 ≥ 1 |
| `dada2.quality_trim_score` | `25` | Illumina：自动估计剪切位置时，平滑后的第 25 百分位 Phred 质量需达到的阈值 | 整数 `0–93` |

PacBio 最短长度不得大于最长长度。`quality_trim_score` 只用于 Illumina 的质量分布剪切位置估计，不是对全部平台逐条 read 的平均质量过滤阈值，也不作用于 Ion Torrent 或 PacBio 分支。

**Ion Torrent 自动 EE**

两种方法默认对当前数据集去引物后的全部 reads 统计长度，按 read 数加权计算精确中位数 `L`，使用 `EE = L × 10^(-19/10)`。例如 `L=407 bp` 时，阈值约为 `5.123826`。统计发生在方法专用的额外前端剪切和质量过滤之前；不会先过滤再反复估计。

同一数据集的 reads 使用一个固定 EE 阈值；这不等于逐 read 按 Q19 过滤，也不会在遇到 Q19 碱基时截断。dada2 仍保留默认的 `trunc-q=2`。该自动值是一项过滤设置，不代表测得的真实错误数。

结果目录中的 `ion_quality-vsearch.json` / `ion_quality-dada2.json` 记录输入数量、中位长度、自动值及实际采用值。要固定使用 EE=2，可设置 `{"vsearch": {"ion_maxee": 2}, "dada2": {"ion_maxee": 2}}`；省略字段或设为 `null` 恢复自动计算，`0` 表示严格的 EE=0。

**ONT 的 `vsearch` 分支：`ont`**

程序按数据集的读长分布估计中心长度：优先使用明显主峰，主峰不明显时使用中位数。读长窗口通常为中心的 `±length_tolerance`，下限至少为 `length_floor`。

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `ont.quality` | `10` | chopper 过滤的最低 read 质量分数 | 整数 `0–93` |
| `ont.length_tolerance` | `0.15` | 自动读长窗口相对于中心的半宽比例；`0.15` 为上下各 15% | 数值 `0–1` |
| `ont.length_floor` | `200` | 自动读长窗口的硬性下限，单位 bp | 整数 ≥ 1 |
| `ont.cluster_identity` | `0.97` | 校正后序列的最终聚类相似度比例 | 数值 `0.01–1` |
| `ont.map_identity` | `0.90` | 过滤后的 reads 回贴到最终代表序列的相似度比例 | 数值 `0.01–1` |

**合并、分类与建树：`taxa`**

| JSON 字段 | 默认值 | 作用 | 取值范围 |
|---|---:|---|---|
| `taxa.confidence` | `0.7` | 分类器置信阈值；命令行对应 `--confidence` | 数值 `0–1` |
| `taxa.singleV` | `false` | 使用单一区域 MAFFT/FastTree 从头建树；命令行对应 `--singleV` | 布尔值 |
| `taxa.notree` | `false` | 跳过建树及树相关过滤；命令行对应 `--notree` / `--no-tree` | 布尔值 |

`singleV` 和 `notree` 不能同时为 `true`；两者均为 `false` 时使用 SEPP 插入建树。分类器种类仍通过命令行 `--classifier greengenes` 或 `--classifier silva` 选择。

**配置示例**

将以下内容保存为 `settings.json`：

```json
{
  "primer": {
    "skip_unknown": true
  },
  "vsearch": {
    "min_frequency": 3
  },
  "taxa": {
    "confidence": 0.8
  }
}
```

例如，本地清单 `local_metadata.csv` 为：

```csv
datasets,path,platform
bee_local,local_data/bees,ILLUMINA
```

其中 `local_data/bees` 相对于 CSV 所在目录。数据集名、路径和平台三个列名均须显式指定。数据与清单准备好后，从项目工作目录运行：

```bash
Meta2Data AmpliconPIP \
  --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path \
  --local-platform-colNAME platform \
  --vsearch \
  --parameter settings.json \
  -t 8

Meta2Data AmpliconTAXA \
  --vsearch \
  --parameter settings.json \
  --notree \
  -t 8
```

PIP 使用引物及 `vsearch` 设置，结果写入 `results/pip/`；TAXA 读取这些结果，使用置信阈值 `0.8`，并按命令行要求跳过建树，输出到 `results/taxa/`。

所有字段的完整默认配置见 [parameters.default.json](parameters.default.json)。

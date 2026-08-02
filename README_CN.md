# Tn5cope

[English](README.md) | **中文**

> Clear signals, careful calls.

**Tn5cope** 是一个面向 Tail-PCR/Sanger 序列的 Tn5 插入位点分析工具。它会尽量保守地
定位 Tn5–genome junction，注释插入点所在基因或基因间区，并输出候选基因、位置簇、
功能簇和基因组圈图。

仓库自带一套完全人工合成的示例输入，可以先用它确认安装正确，再换成自己的文件。

> [!IMPORTANT]
> 本公开仓库不包含原始实验序列、私有参考基因组文件或表格形式的分析结果。
> `examples/` 中的输入全部是人工合成的教学数据，不能用于任何生物学结论。下图仅用于
> 展示预期输出，是不含原始序列的衍生可视化结果。

<p align="center">
  <img src="assets/representative_genome_circle_plot.png" width="760" alt="Tn5cope 生成的代表性基因组圈图">
</p>
<p align="center"><em>Tn5cope 完成一次分析后生成的代表性结果。</em></p>

## 它能做什么

- 在 query 的正向和反向互补方向搜索可靠的基因组定位。
- 识别 Tn5 mosaic end，并以 Tn5–genome junction 相邻碱基作为插入坐标。
- 将插入点区分为 `CDS`、`gene_non_CDS` 或 `intergenic`。
- 为基因间插入提供方向、距离、operon 和 Tn5 构建体证据驱动的候选排序。
- 将失败、低置信和多重命中单独保留，不强行给出一个“最好看”的答案。
- 输出可追溯的 TSV、Excel、位置簇、功能簇、热图和基因组圈图。

Tn5cope 适合单个突变体的 Tail-PCR/Sanger 插入位点定位，**不是**用于 Tn-seq 丰度、
基因必需性、显著性富集或表型因果证明的软件。

## 第一次运行：先用合成示例

### 1. 准备 Python

需要 Python 3.9 或更高版本。在终端中检查：

```bash
python3 --version
```

如果系统显示 `Python 3.9` 或更高版本，就可以继续。Windows 上的命令可能是
`python --version`。

### 2. 进入项目目录并创建独立环境

在 macOS/Linux 终端中：

```bash
cd /path/to/Tn5cope
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

在 Windows PowerShell 中：

```powershell
cd C:\path\to\Tn5cope
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

命令行最前面出现 `(.venv)`，通常表示环境已经激活。

### 3. 运行教学示例

下面的文件全部位于 [`examples/demo`](examples/demo)，没有任何真实实验数据：

```bash
tn5cope \
  --query-csv examples/demo/query.csv \
  --genome-fasta examples/demo/reference.fna \
  --gtf examples/demo/reference.gtf \
  --operon-tsv examples/demo/operons.tsv \
  --aligner internal \
  --output-dir demo_output \
  --output-prefix demo
```

Windows PowerShell 用户可以把上面的命令写在同一行运行。完成后，打开新生成的
`demo_output/`，重点查看：

- `demo_main_results.tsv`：每条 query 的定位与注释结果；
- `demo_failed_or_ambiguous_hits.tsv`：没有可靠定位的记录；
- `demo.xlsx`：适合用 Excel 浏览的汇总工作簿；
- `demo_genome_circle_plot.png`：位置分布圈图。

如果看到这些文件，说明安装和基本流程正常。`demo_output/` 已被 `.gitignore` 排除，
不会被上传。

## 准备自己的输入

正式运行需要三个文件，另有一个可选文件：

| 文件 | 是否必需 | 作用 |
|---|---|---|
| query CSV | 必需 | Tail-PCR/Sanger 得到的待定位序列 |
| reference FASTA | 必需 | 与实验菌株匹配的参考基因组 |
| reference GTF | 必需 | 与 FASTA 同一版本的 gene/CDS 注释 |
| operon TSV | 可选 | 已知或外部预测的 operon 关系 |

最安全的做法是把真实文件放在仓库外，例如单独的 `private_input/` 文件夹，然后在命令中
填写它们的完整路径。不要用真实文件覆盖 `examples/`。

### 1. Query CSV

推荐使用单列 CSV，第一行写 `sequence`，之后每行一条序列：

```csv
sequence
CTGTCTCTTATACACATCTACGTTGCAAGTCCGATGCTAACGGTTCAGATC
GGAACCTTAGCGTTCGATACCGTTAAGCTGACCTAGTCCGAT
```

要求与注意事项：

- 第一列是序列；当前 CLI 不读取样本名列，而是依次生成 `Q001`、`Q002`……
- 建议保留表头；没有表头也可以读取。
- 大小写均可，空白会被清理；建议只使用 `A/C/G/T/N`。
- 序列可以包含完整或部分 Tn5 mosaic end，也可以只包含基因组侧翼。
- 如果没有检测到 Tn5 motif，程序会把输入序列的**第一个碱基**视为靠近插入端的一侧。
  因此必须先确认导出的 Sanger 序列方向符合自己的实验设计。
- 去除 Tn5 后的可比对侧翼不能短于 20 bp；真实数据通常应尽量提供更长的高质量序列。

[`examples/demo/query.csv`](examples/demo/query.csv) 是可以直接打开查看的人工模板。

### 2. Reference FASTA

FASTA 的第一行以 `>` 开头，后面是参考序列：

```fasta
>chromosome_name
ACGTTGCAAGTCCGATGCTAACGGTTCAGATCGTACGACT...
```

当前主流程按**单 contig/单染色体**设计。多 contig 文件会报错，不应把多个 contig
简单拼接。FASTA 头部第一个空格前的名称，例如 `chromosome_name`，必须与 GTF 第一列
完全一致。

### 3. Reference GTF

GTF 是以 Tab 分隔的九列表格。Tn5cope 至少需要 `gene` 和 `CDS` 记录：

```text
chromosome_name  source  gene  101  900  .  +  .  gene_id "GENE001"; gene "demoA"; locus_tag "LOC_0001";
chromosome_name  source  CDS   101  900  .  +  0  gene_id "GENE001"; gene "demoA"; locus_tag "LOC_0001"; product "example protein";
```

上面为了阅读使用空格展示；实际文件的九列必须使用 Tab。应特别检查：

- GTF 与 FASTA 来自同一个参考版本；
- GTF 第一列和 FASTA 序列名称一致；
- 坐标是合法整数，链方向为 `+` 或 `-`；
- 属性中尽量包含 `gene_id`、`locus_tag`、`gene` 和 `product`；
- Ontology/GO 信息不是定位必需项，但会影响功能聚类的可解释性。

### 4. Operon TSV（可选）

如果有可靠的 operon 注释，可以提供 Tab 分隔文件：

```tsv
operon_id	locus_tag	gene_name
OP001	LOC_0001	demoA
OP001	LOC_0002	demoB
```

必须包含 `operon_id`，并至少包含一个基因标识列：`locus_tag`、`gene_id`、
`gene_code`、`gene_name` 或 `gene`。没有可靠注释时直接省略 `--operon-tsv`，程序不会
因为两个基因相邻就擅自宣布它们属于同一 operon。

## 运行自己的数据

将路径替换为自己的真实文件：

```bash
tn5cope \
  --query-csv /path/to/private_input/query.csv \
  --genome-fasta /path/to/private_input/reference.fna \
  --gtf /path/to/private_input/reference.gtf \
  --tn5-structure-profile unknown \
  --output-dir /path/to/private_output \
  --output-prefix tn5cope_run
```

`--aligner auto` 会优先寻找 minimap2、BLASTn 或 BWA；都不可用时使用内置的可复现
fallback。教学示例可以使用 `--aligner internal`，正式数据更推荐安装至少一个成熟
比对器，例如 minimap2。macOS/Homebrew 用户可以运行：

```bash
brew install minimap2 blast bwa
```

这些比对器是系统程序，不会随 `pip install` 自动安装。也可以显式选择
`--aligner minimap2`、`--aligner blastn`、`--aligner bwa` 或
`--aligner internal`。

如果实际 Tn5 构建体结构不确定，请保留
`--tn5-structure-profile unknown`。只有在构建体图谱明确证明存在外向启动子或转录
终止子时，才选择其他 profile；该选项会改变基因间候选排序。

## 如何理解主要输出

- `*_main_results.tsv`：最终插入坐标、定位置信状态和直接基因注释。
- `*_raw_matches.tsv`：候选命中依据；这里的 “raw” 指未汇总的**比对记录**，不是原始
  Sanger 序列。
- `*_failed_or_ambiguous_hits.tsv`：无命中、低置信或歧义记录。
- `*_expanded_gene_events.tsv`：正式 Expanded candidate 的 query–gene 事件。
- `*_expanded_candidate_exclusions.tsv`：未进入 Expanded 范围的事件及原因。
- `*_expanded_functional_cluster_*.tsv`：功能簇成员和特征摘要。
- `*_position_cluster_*.tsv`：全部高置信位置的 circular DBSCAN 结果。
- `*_expanded_position_cluster_*.tsv`：Expanded candidate 范围的位置簇。
- `*_genome_circle_plot.png/.pdf/.svg`：基因组圈图。
- `*.xlsx`：适合人工检查的多 sheet 工作簿。

`Insertion site` 是与 Tn5–genome junction 直接相邻的 1-based 参考碱基；
`Junction_boundary_0based` 是两个参考碱基之间的 0-based 边界。`Feature_relation`
描述插入点直接落在何处；`Candidate ...` 列是基于方向、距离和可用注释作出的初步推断。
两者不要混为一谈。

功能簇表示现有注释的相似性，位置簇表示给定 DBSCAN 参数下的操作性分组；它们都不是
富集显著性或因果证明。候选基因仍需 genome-specific PCR/Sanger、互补、定点突变或
表达实验验证。

## 可选的本地网页

不习惯命令行时，可以安装本地网页依赖：

```bash
python -m pip install -e ".[web]"
tn5cope-web
```

浏览器打开 <http://127.0.0.1:8000>。网页仍在自己的电脑上运行，上传文件会暂存在
仓库下的 `runs/`；该目录不会被 Git 上传，但不会自动加密，使用后请按实验室数据管理
规范处理。

## 常见问题

### `tn5cope: command not found`

通常是虚拟环境没有激活。重新进入项目目录并运行 `source .venv/bin/activate`；Windows
PowerShell 使用 `.\.venv\Scripts\Activate.ps1`。

### 大部分序列都是 `no_hit`

先检查参考菌株/版本、query 的插入端方向、侧翼长度和 Sanger 质量。然后查看
`*_raw_matches.tsv` 与 `*_failed_or_ambiguous_hits.tsv`，不要仅靠放宽阈值强行定位。

### GTF 中明明有基因，结果却没有注释

最常见原因是 FASTA 标题与 GTF 第一列不一致，或者 GTF 中缺少 `gene`/`CDS` 记录。

### 能否分析多 contig 组装？

当前正式流程不能。不要把多个 contig 直接拼接成一条假染色体，否则会制造错误坐标和
错误的环状距离。

## 数据隐私

真实输入和输出应始终放在克隆仓库之外。项目自带的 `.gitignore` 会阻止常见序列、参考
文件和结果格式被普通 Git 操作加入版本库，但它不是加密工具，也无法保护手动强制添加的
文件。只有 `examples/` 下的人工合成模板适合纳入版本控制。

## License

Tn5cope 以 [MIT License](LICENSE) 开源。你可以使用、修改和分发本项目，但必须保留原始
版权与许可声明。软件按现状提供，不附带任何担保。

# SIF / Apptainer 使用说明

使用包含当前命令行接口的新版 SIF，直接运行镜像内的 Meta2Data 和依赖。下载命令与 `m2d` 快捷命令配置见 [README](../README.md#sif--apptainer)。

在工作目录中准备 `metadata.csv`、`local_metadata.csv` 和 `local_data/`，然后运行命令。结果与参考资源固定保存到：

```text
工作目录/results/pip/     # 在线与本地测序数据处理结果
工作目录/results/taxa/    # 合并、分类注释及建树结果
工作目录/results/db/      # 自动准备的参考资源
```

TAXA 可用 `-i` 指定已有 PIP 结果，输出仍固定为当前工作目录下的 `results/taxa/`。不同分析项目应切换工作目录。MetaDL 的元数据输出仍通过 `-o` 单独指定。

**直接运行**

```bash
apptainer exec --cleanenv "$M2D_SIF" Meta2Data AmpliconPIP \
  --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path \
  --local-platform-colNAME platform \
  --dada2 -t 8

apptainer exec --cleanenv "$M2D_SIF" Meta2Data AmpliconTAXA \
  --dada2 --classifier greengenes --notree -t 8
```

示例使用 README 中配置的绝对路径变量 `M2D_SIF`，切换分析目录后仍指向同一个镜像。

**使用 `m2d` 快捷命令**

```bash
m2d AmpliconPIP \
  --public-m metadata.csv \
  --public-bioproject-colNAME Bioproject \
  --public-sra-colNAME Run \
  --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path \
  --local-platform-colNAME platform \
  --dada2 -t 8

m2d AmpliconTAXA --dada2 --classifier greengenes -t 8
```

镜像包含程序及依赖，无需挂载源码或个人软件环境。QIIME 配置、Numba 和绘图缓存使用可写的临时目录；设置 `TMPDIR` 时，该目录须在容器内可访问。

本地 CSV 的相对路径以 CSV 所在目录为基准。数据集名、路径和平台三列均须通过对应参数显式指定，没有默认列名；即使列名是 `datasets`、`path`、`platform` 也要填写。需要读取引物列时，再添加 `--local-primer-f-colNAME` 和可选的 `--local-primer-r-colNAME`。

**内置测试**

```bash
m2d AmpliconPIP --test --vsearch -t 8
m2d AmpliconTAXA --vsearch --notree -t 8
```

内置清单随镜像提供；运行测试会下载真实 Run。扩展清单的获取方式见 [README 测试集](../README.md#test-datasets)。

## 在 GitHub 网页上构建并发布新版 SIF

`Build SIF` 工作流在 GitHub Actions 中先验证原生 Conda 环境，再从 Ubuntu 24.04 安装新的 QIIME 2 2026.7 环境、构建 SIF，并测试镜像内的实际程序。全部通过后才发布到 GitHub Container Registry（GHCR）。工作流只由 `workflow_dispatch` 手动触发，不会因普通提交或 PR 自动运行。

兼容的 GG2 2024.09 和 SILVA 138.1 模型已通过核验并发布。原生环境通过全部 317 项单元测试（零失败、错误或跳过）、14 个 PIP 场景和 5 个 TAXA 场景。最终 SIF 通过 2 个 Illumina PE PIP 场景（DADA2、VSEARCH）和 5 个 TAXA 场景。源码版本、范围及证据见 [最终验证记录](validation/2026-11-19/release-gates.md)。

### 首次配置

确认仓库已启用 Actions，默认分支存在支持 `workflow_dispatch` 的 `Build SIF` 工作流。本仓库原有工作流已满足这一条件。升级代码和新版工作流从下方 **Use workflow from** 选择的候选分支运行，完成原生环境与最终 SIF 验证后再合并。在现有 `meta2data` 包的 **Package settings → Manage Actions access** 中，允许 `LinyangSun/Meta2Data` 仓库写入。工作流使用自动提供的 `GITHUB_TOKEN`，权限为 `contents: read` 和 `packages: write`，无需创建个人 Token。

### 网页操作

1. 打开仓库的 **Actions → Build SIF → Run workflow**。
2. 在 **Use workflow from** 中选择包含待发布代码的分支。工作流打包触发时该分支的提交，不会自动合并其他 PR。
3. 填写 `output_tag`：使用尚未占用的新标签，不能使用 `latest`。工作流检查目标标签是否存在，已存在则停止，不覆盖旧标签。新构建不再需要 `base_tag`。
4. 点击 **Run workflow**。原生环境验证通过后才开始构建；镜像验证通过后才上传。运行日志和 validation artifacts 保存具体检查结果。
5. 在运行摘要中查看源码提交、新 SIF 的 SHA-256，以及可直接复制的下载命令。

正式版本 `v1.1` 已发布到 [GitHub Packages](https://github.com/LinyangSun/Meta2Data/pkgs/container/meta2data)。在服务器下载为独立文件：

```bash
apptainer pull Meta2Data-v1.1.sif \
  oras://ghcr.io/linyangsun/meta2data:v1.1

export M2D_SIF="$(pwd)/Meta2Data-v1.1.sif"
```

更新 `M2D_SIF` 指向新文件后，`m2d` 会使用该镜像。

### 构建范围

`container/Meta2Data.def` 使用 `Bootstrap: docker` 和 Ubuntu 24.04，通过带校验和的 Miniforge 安装器创建 `/opt/conda/envs/meta2data`。`env.yml` 只列出实际调用的 8 个 QIIME 插件、核心程序及其他直接依赖，QIIME 2 固定为 2026.7.0、scikit-learn 固定为 1.7.1。`environment-constraints-linux-64.txt` 保留官方发布环境中实际需要的依赖版本和构建，不重新加入未使用的插件。Conda 用户与 SIF 使用同一个 `environment-linux-64.lock`，锁定这些包的下载地址和 SHA256；镜像中的实际安装清单另存为 `/opt/meta2data-current/environment-installed-linux-64.txt`。

原生安装、CI 与 SIF 构建统一调用 `scripts/create_conda_environment.py`。脚本先核对 `env.yml` 和依赖约束文件的 SHA256，再逐项核对锁文件中的完整包集合、版本和构建。随后通过 `conda create --override-channels --no-default-packages` 安装，不加入用户配置的额外 channel 或默认包。这样可避开 Conda 26.7–26.9 对 `nodefaults` 的已知问题；不会修改用户的 Conda 配置。

环境不安装 `q2-fondue`、Conda `sra-tools` 或 `ossuuid`。同一安装脚本会将官方 NCBI SRA Toolkit 3.4.1 的固定校验和归档安装到 `/opt/conda/envs/meta2data/opt/sratoolkit.3.4.1-alma_linux64/`，在环境内提供 `prefetch` 和 `fasterq-dump`。该目录的 `meta2data-install.json` 保存来源、校验和与工具版本，CI 同时保存这份记录。这样避免了已确认的 SRA 可执行文件重复归属，以及 `ossuuid` 和 `libuuid` 的文件冲突。

SEPP 参考仍安装在 `/opt/Meta2Data/reference/`。三个分类器使用各自固定版本的资源文件和兼容性检查，按需准备到工作目录的 `results/db/`。

构建工具固定为 Apptainer 1.5.5，构建时执行 `%test`。测试顺序为：

1. 原生 Conda 环境：版本、插件、工具、单元测试及小型数据的完整平台分支验证。
2. SIF 构建中的 `%test`：确认镜像内实际安装的版本和命令入口。
3. 最终 SIF：运行小型 Illumina PE 的 DADA2/VSEARCH、三个分类器、SEPP 和从头建树检查。

任何一步失败都会阻止发布。镜像检查直接运行镜像内的代码和依赖，不挂载宿主 Conda 环境。兼容性检查不评估数据库注释的生物学准确性。

### 云端磁盘空间

默认使用 `ubuntu-24.04` GitHub-hosted runner。工作流仅在默认的临时 runner 上移除本任务不使用的预装 Android、.NET 和 GHC 目录来释放空间，不操作你的服务器或个人文件。

新环境、解压的容器文件系统、参考资源、测试结果及最终 SIF 会同时占用空间。构建前要求至少有 30 GiB 可用空间；这是运行门槛，不是最终占用的保证。空间不足时，可将仓库 **Settings → Secrets and variables → Actions → Variables** 中的 `SIF_RUNNER` 设为该仓库有权使用的、更大磁盘的 Ubuntu x64 GitHub-hosted runner 标签。此变量不会创建机器或开通付费服务；使用 larger runner 可能产生费用。

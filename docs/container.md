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
apptainer exec --cleanenv Meta2Data.sif Meta2Data AmpliconPIP \
  --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path \
  --local-platform-colNAME platform \
  --dada2 -t 8

apptainer exec --cleanenv Meta2Data.sif Meta2Data AmpliconTAXA \
  --dada2 --classifier greengenes --notree -t 8
```

示例假设镜像在当前目录；位于其他位置时，将 `Meta2Data.sif` 换成镜像路径。

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

`Build SIF` 工作流在 GitHub Actions 云端拉取旧 SIF、打包所选版本的程序代码，然后直接发布到 GitHub Container Registry（GHCR）。不需要本地打包、上传，也不需要 Codespaces。工作流只由 `workflow_dispatch` 手动触发，不会因普通提交或 PR 自动运行。

### 首次配置

将 `.github/workflows/build-sif.yml` 合并到默认分支，并确保仓库已启用 Actions。在现有 `meta2data` 包的 **Package settings → Manage Actions access** 中，允许 `LinyangSun/Meta2Data` 仓库写入。工作流使用自动提供的 `GITHUB_TOKEN`，权限为 `contents: read` 和 `packages: write`，无需创建个人 Token。

### 网页操作

1. 打开仓库的 **Actions → Build SIF → Run workflow**。
2. 在 **Use workflow from** 中选择包含待发布代码的分支。工作流打包触发时该分支的提交，不会自动合并其他 PR。
3. 填写 `base_tag`：GHCR 中已有、可正常使用的 Meta2Data SIF 标签。可在 [GitHub Packages](https://github.com/LinyangSun/Meta2Data/pkgs/container/meta2data) 中查看。
4. 填写 `output_tag`：尚未使用的新标签，例如 `update-20261028`。不能与基础标签相同，也不能使用 `latest`。工作流检查目标标签是否存在，已存在则停止，不覆盖旧标签。
5. 点击 **Run workflow**。完成后在运行摘要中查看源码提交、基础和新 SIF 的 SHA-256，以及可直接复制的下载命令。

例如，发布标签为 `update-20261028` 后，在服务器下载为新文件：

```bash
apptainer pull Meta2Data-update-20261028.sif \
  oras://ghcr.io/linyangsun/meta2data:update-20261028
```

上面的标签仅为命名示例，不表示该镜像已经发布。更新 `M2D_SIF` 指向新文件后，`m2d` 才会使用新版镜像。

### 构建范围

工作流沿用 `container/Meta2Data.def` 的 `Bootstrap: localimage`，保留基础镜像中的 QIIME 2 环境、其他软件依赖和 SEPP 参考资源，只更新定义文件列出的程序目录和文件。`env.yml` 会复制到镜像中，但不会执行 Conda 环境更新；增加或升级依赖需要另行更新基础镜像或构建定义。

构建工具固定为 Apptainer 1.5.5；构建命令使用 `--notest`，不执行定义文件的 `%test`、单元测试或分析测试。构建和上传成功只表示生成了新镜像，不表示分析流程已通过运行验证。该工作流不修改分析代码、检查点机制或旧镜像标签。

### 云端磁盘空间

默认使用 `ubuntu-24.04` GitHub-hosted runner。工作流仅在默认的临时 runner 上移除本任务不使用的预装 Android、.NET 和 GHC 目录来释放空间，不操作你的服务器或个人文件。

下载基础 SIF 后，工作流以基础镜像大小的四倍加 2 GiB 估算后续构建所需空间。该估算不等于实际解压大小；通过检查也不能保证大镜像一定能构建成功。镜像下载或解压时空间不足，可将仓库 **Settings → Secrets and variables → Actions → Variables** 中的 `SIF_RUNNER` 设为该仓库有权使用的、更大磁盘的 Ubuntu x64 GitHub-hosted runner 标签，然后重新运行。此变量不会创建机器或自动开通付费服务；使用 larger runner 可能产生费用。

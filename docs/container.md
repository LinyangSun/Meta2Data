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

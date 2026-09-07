# scPhenomics 分析流水线
[![English](https://img.shields.io/badge/English-blue?style=flat-square)](./../../README.md)[![中文](https://img.shields.io/badge/中文-red?style=flat-square)](README_zh-cn.md)
**一套用于多通道荧光显微镜的通用单细胞表型组图像分析流水线**
## 概述
本代码仓库是 [CCRMicroscopyCore/kangs](https://github.com/CCRMicroscopyCore/kangs)（Kangs 等人，*Nature Metabolism*, 2025）的复刻分支与通用扩展版本。
⚠️ **注意**：该复刻分支面向通用成像表型组分析开发，**不会定期与上游仓库同步**。
如果需要论文对应的原始实现，请查阅[上游仓库](https://github.com/CCRMicroscopyCore/kangs)。
本流水线源自 Kang S. 等人发表于《Nature Metabolism》2025年的工作，提供肝小叶定量与单细胞表型组分析的完整工作流。
## 环境配置
### 1. 安装 Python 3.13
确保系统已安装 [Python 3.13](https://www.python.org/)。
### 2. 安装 Poetry（推荐）
[Poetry](https://python-poetry.org/docs/) 是 Python 的依赖管理与打包工具，可以声明项目依赖，便捷完成依赖的安装/更新，同时提供锁定文件保证可复现安装。
如果你更倾向使用 `uv`、`pip` 或 `conda` 等其他工具，请参考 [`pyproject.toml`](pyproject.toml) 文件查看所需库的版本要求。
### 3. 安装 VSCode（可选但推荐）
建议安装 [VSCode](https://code.visualstudio.com/) 并搭配 **Python 语言支持扩展**。针对数据科学工作流，可查阅 VSCode 数据科学文档获取更多扩展。
此处不提供详细安装步骤，请查阅各工具官方文档。
## 克隆代码仓库
建议先 Fork 到自己的 GitHub 账号，也可以直接克隆：
```bash
git clone https://github.com/Mcatstar/scPhenomics-pipeline.git
cd scPhenomics-pipeline
```
## 创建虚拟环境并安装依赖
### 1. 配置 Poetry，在项目根目录生成 `.venv`
```bash
poetry config virtualenvs.in-project true
```
### 2. 查看可用 Python 版本
```bash
poetry python list
```
### 3. 指定 Python 版本并创建虚拟环境
```bash
poetry env use /path/to/python3.13
```
### 4. 安装依赖
```bash
poetry install --no-root
```
⚠️ 本项目仅为脚本流水线，**并非可安装的 Python 软件包**。
`--no-root` 参数用于阻止 Poetry 将当前项目作为软件包进行安装。
### 在 VSCode 打开项目
```bash
code .
```
## 使用说明
1. 使用 `notebooks/` 目录下的 Jupyter 笔记本**验证数据分析工作流**。
2. 工作流验证通过后，在 `src/` 目录编写正式生产脚本。
### 图像分割与定量分析
图像分割与定量分析基于 **Python 3.13** 实现。
- 分割使用定制的 [Cellpose](https://github.com/MouseLand/cellpose) 模型。
- 建议在开启 GPU 的高性能计算集群运行分割任务，但并非强制要求。
- 图像可视化与校验推荐使用 [napari](https://github.com/napari/napari)（多维图像查看器），也可以使用 FIJI 等其他工具。
#### `lobule_analysis_v2.ipynb`
该笔记本执行以下分析步骤：
- **加载图像**
  定义肌动蛋白、线粒体、脂质各个成像通道。
- **细胞分割**
  结合肌动蛋白与线粒体通道完成分割。
  → 支持手动修改细胞标签。
- **线粒体分割**
- **脂质分割**
- **生成重叠映射图**
  将线粒体与脂质的分割结果叠加。
- **生成欧氏距离变换（EDT）映射图**
  - 基于人工标记的中央静脉
  - 基于人工标记的门管区静脉
- **生成基于细胞器的距离变换映射图**
  - 基于线粒体
  - 基于脂质
- **定量统计**
  针对细胞、线粒体、脂质、细胞器重叠区域统计：
  - 通道信号强度
  - 面积
  - 质心坐标
  - 几何形态参数
### 数据分析
数据分析脚本基于 **R (v4.2)** 编写。
#### `PLIN5_analysis.R`
- 读取定量结果数据
- 计算细胞与细胞器相对中央静脉‑门管区静脉距离，归一化至中央‑门管最大距离
- 根据相对距离将数据划分为 **12个区间（R1‑R12）**：
  - **R1** → 靠近门管区静脉(PV)
  - **R12** → 靠近中央静脉(CV)
## 许可证
本项目基于 **MIT 许可证**开源。
原始版权归原作者所有，完整协议请查看 [LICENSE](https://license/) 文件。
## 引用
如果你使用该流水线，请引用：
- **原始方法与预训练 Cellpose 权重**：
  Kang S. et al., *Nature Metabolism*, 2025，以及对应的[上游仓库](https://github.com/CCRMicroscopyCore/kangs)。
## 贡献
欢迎提交 Issue 与 Pull Request。
## 联系
如有疑问或建议，请在本仓库提交 Issue。
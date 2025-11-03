import pandas as pd
import os
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

# ====== 通用配置 ======
base_dir = "/public/home/Chaochen/cellphonedb"
cpdb_file_path = f"{base_dir}/v5.0.0/cellphonedb.zip"
mapping_path = f"{base_dir}/mart_export.txt"

# ====== 定义各组的路径 ======
datasets = {
    "ND": {
        "counts": f"{base_dir}/ND_counts.txt",
        "meta": f"{base_dir}/ND_meta.tsv",
        "out": f"{base_dir}/result_cpdb_DEG_ND"
    },
    "HFD": {
        "counts": f"{base_dir}/HFD_counts.txt",
        "meta": f"{base_dir}/HFD_meta.tsv",
        "out": f"{base_dir}/result_cpdb_DEG_HFD"
    }
}

# ====== 载入基因名映射表 ======
print("🔁 正在加载本地基因映射表...")
mapping = pd.read_csv(mapping_path, sep="\t")
mapping = mapping.dropna(subset=["Gene name", "Human gene name"])
gene_map = dict(zip(mapping["Gene name"], mapping["Human gene name"]))

# ====== 分别处理每组数据 ======
for group, paths in datasets.items():
    print(f"\n=== 🚀 开始处理: {group} ===")

    # --- Step 1: counts 映射 ---
    print(f"📂 加载 counts: {paths['counts']}")
    counts = pd.read_csv(paths["counts"], sep="\t", index_col=0, encoding='utf-8', engine='python')
    counts_human = counts[counts.index.isin(gene_map)]
    counts_human.index = counts_human.index.map(gene_map)
    counts_human = counts_human[~counts_human.index.duplicated()]
    counts_human_file = paths["counts"].replace(".txt", "_human.txt")
    counts_human.to_csv(counts_human_file, sep="\t")
    print(f"✅ 映射完成，保留基因数: {counts_human.shape[0]}")



import pandas as pd
import scanpy as sc
from anndata import AnnData

# ===== 输入路径 =====
nd_counts_path = "/public/home/Chaochen/cellphonedb/ND_counts_human.txt"
nd_meta_path = "/public/home/Chaochen/cellphonedb/ND_meta.txt"
hfd_counts_path = "/public/home/Chaochen/cellphonedb/HFD_counts_human.txt"
hfd_meta_path = "/public/home/Chaochen/cellphonedb/HFD_meta.txt"

# ===== 输出路径 =====
nd_h5ad_path = nd_counts_path.replace(".txt", ".h5ad")
nd_meta_out_path = nd_meta_path.replace(".txt", ".tsv")
hfd_h5ad_path = hfd_counts_path.replace(".txt", ".h5ad")
hfd_meta_out_path = hfd_meta_path.replace(".txt", ".tsv")

# ===== 函数定义：读取 counts 和 meta，生成 h5ad 和 tsv =====
def process_group(counts_path, meta_path, h5ad_path, meta_out_path, group_name):
    print(f"\n🔁 正在读取 {group_name} counts...")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)

    print(f"💾 {group_name} 保存为 .h5ad 格式...")
    adata = AnnData(X=counts.T)  # Cell × Gene
    adata.var_names = counts.index
    adata.obs_names = counts.columns
    adata.write(h5ad_path)
    print(f"✅ 已保存: {h5ad_path}")

    print(f"📋 正在转换 {group_name} meta 为 .tsv...")
    meta = pd.read_csv(meta_path, sep="\t")
    meta.to_csv(meta_out_path, sep="\t", index=False)
    print(f"✅ meta 保存为: {meta_out_path}")

# ===== 执行处理 =====
process_group(nd_counts_path, nd_meta_path, nd_h5ad_path, nd_meta_out_path, group_name="ND")
process_group(hfd_counts_path, hfd_meta_path, hfd_h5ad_path, hfd_meta_out_path, group_name="HFD")

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# === 文件路径定义 ===
base_dir = "/public/home/Chaochen/cellphonedb"
cpdb_file_path = f"{base_dir}/v5.0.0/cellphonedb.zip"

datasets = {
    "ND": {
        "counts": f"{base_dir}/ND_counts.h5ad",
        "meta": f"{base_dir}/ND_meta.tsv",
        "out": f"{base_dir}/result/statistical_analysis_ND"
    },
    "HFD": {
        "counts": f"{base_dir}/HFD_counts.h5ad",
        "meta": f"{base_dir}/HFD_meta.tsv",
        "out": f"{base_dir}/result/statistical_analysis_HFD"
    }
}

# === 运行 CellPhoneDB 分析函数 ===
def run_cpdb(meta_path, counts_h5ad, out_path, group_name):
    print(f"\n🚀 正在运行 CellPhoneDB 分析: {group_name}")
    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path=cpdb_file_path,
        meta_file_path=meta_path,
        counts_file_path=counts_h5ad,
        counts_data='hgnc_symbol',
        score_interactions=True,
        iterations=1000,
        threshold=0.1,
        threads=5,
        debug_seed=42,
        result_precision=3,
        pvalue=0.05,
        subsampling=False,
        subsampling_log=False,
        subsampling_num_pc=100,
        subsampling_num_cells=1000,
        separator='|',
        debug=False,
        output_path=out_path,
        output_suffix=None
    )
    print(f"✅ {group_name} 分析完成，结果保存在：{out_path}")
    return cpdb_results

# === 分别运行 ND 和 HFD ===
for group, paths in datasets.items():
    run_cpdb(paths["meta"], paths["counts"], paths["out"], group)

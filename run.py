import argparse
import scanpy as sc
import numpy as np
from pathlib import Path
import pandas as pd
import json



# Read dataset
parser = argparse.ArgumentParser()
parser.add_argument("--output_dir", "-o", required=True)
parser.add_argument("--name", "-n", required=True)
parser.add_argument("--data.h5ad", dest="data_h5ad", required=True)
parser.add_argument(
    "--data.clusters_truth_num",
    dest="clusters_truth_num_path",
    required=True,
)
parser.add_argument("--d_cluster", type=int, default=0)
parser.add_argument("--seed", type=int, default=42)

args, _ = parser.parse_known_args()

adata = sc.read_h5ad(args.data_h5ad)

with open(args.clusters_truth_num_path, "rb") as f:
    n_clusters = int(f.read().strip()) + args.d_cluster

# Generate random labels
random_range = np.random.default_rng(args.seed)
labels = random_range.integers(0, n_clusters, size=adata.n_obs)

# Output directory
out = Path(args.output_dir)
out.mkdir(parents=True, exist_ok=True)

# Write compatible cluster assignments
clusters = pd.DataFrame({
    "cell_id": adata.obs_names,
    "louvain": labels.astype(str),
    "leiden": labels.astype(str),
})

clusters.to_csv(out/f"{args.name}.clusters.csv", sep="\t", index=False)

# Provide any outputs required by the benchmark
with open(out/f"{args.name}.timings.tsv", "w") as f:
    json.dump({}, f)

with open(out/f"{args.name}.clustering_info.tsv", "w") as f:
    json.dump({}, f)

pd.DataFrame({"cell_id": adata.obs_names}).to_csv(
    out / f"{args.name}.pca.tsv", sep="\t", index=False
)
(out / f"{args.name}.hvgs.tsv").write_text("")
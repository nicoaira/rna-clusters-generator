# RNA Cluster Dataset Generator

This project is a Rust-based RNA cluster dataset generator. Instead of triplets (anchor/positive/negative), it creates clusters with one anchor and multiple positive (modified) structures per cluster. For example, if you choose 50 clusters and 200 structures per cluster, the program will generate 50 anchor structures and for each anchor it will produce 200 “positive” structures by introducing noise.

## Usage Example

You can run the program with:

```bash
cargo run -- \
  --n-clusters 50 \
  --structures-per-cluster 200 \
  --seq-min-len 40 \
  --seq-max-len 600 \
  --seq-len-distribution norm \
  --seq-len-mean 300 \
  --seq-len-sd 100 \
  --n-stem-indels 2 \
  --stem-max-size 14 \
  --stem-min-size 3 \
  --stem-max-n-modifications 1 \
  --n-hloop-indels 2 \
  --hloop-min-size 1 \
  --hloop-max-size 10 \
  --hloop-max-n-modifications 1 \
  --n-iloop-indels 2 \
  --iloop-min-size 2 \
  --iloop-max-size 10 \
  --iloop-max-n-modifications 1 \
  --n-bulge-indels 2 \
  --bulge-min-size 1 \
  --bulge-max-size 8 \
  --bulge-max-n-modifications 1 \
  --n-mloop-indels 2 \
  --mloop-min-size 2 \
  --mloop-max-size 15 \
  --mloop-max-n-modifications 1 \
  --num-workers 12 \
  --batch-size 64 \
  --split \
  --train-fraction 0.97 \
  --val-fraction 0.03
```

## Parameters Explanation

### Cluster Generation
- `--n-clusters`: Number of clusters (each with one anchor).
- `--structures-per-cluster`: Number of positive structures per cluster.

### Sequence Generation and Modification
- (Other parameters control sequence length, stem/loop modifications, etc.)

## Requirements

You can install requirements with

```
chmod +x install_requirements.sh
sudo ./install_requirements.sh
```

### RNAfold
- Install ViennaRNA package from https://www.tbi.univie.ac.at/RNA 
- Ensure that `RNAfold` is accessible in your PATH
- On Ubuntu/Debian, you can install via: `sudo apt install vienna-rna`
- On macOS with Homebrew: `brew install viennarna`
- For other systems, follow the installation guide at https://www.tbi.univie.ac.at/RNA/#download

### Other Requirements
- **ImageMagick**: For merging plots (`montage` command)
- **Rust Toolchain**: Install Rust and Cargo to build and run this project

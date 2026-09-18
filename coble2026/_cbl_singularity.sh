
echo "This script is to download the singularity sif files for the Coble 2026 paper"

echo "Downloading singularity for Panel(b): Reproduce - Scherer et al (2025) ..."
singularity build cbl-scherer-rep.sif docker://ghcr.io/coble-tools/coble:scherer-rep

echo "Downloading singularity for Panel(c): Replicate - Zheng et al (2017) ..."
singularity build cbl-zheng-rep.sif docker://ghcr.io/coble-tools/coble:zheng-rep

echo "Downloading singularity for Panel(d): Generalise - Zheng et al (2017) ..."
singularity build cbl-zheng-gen.sif docker://ghcr.io/coble-tools/coble:zheng-gen


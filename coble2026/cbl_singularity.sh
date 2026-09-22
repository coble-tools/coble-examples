
echo "This script is to download the singularity sif files for the Coble 2026 paper"

echo "Downloading singularity for Panel(b): Reproduce - Scherer et al (2025) ..."
singularity build cbl-scherer-rep.sif docker://ghcr.io/coble-tools/coble:scherer-rep

echo "Downloading singularity for Panel(c): Replicate - Zheng et al (2017) ..."
singularity build cbl-zheng-rep.sif docker://ghcr.io/coble-tools/coble:zheng-rep

echo "Downloading singularity for Panel(d): Replicate-Extend - Zheng et al (2017) ..."
singularity build cbl-zheng-gen.sif docker://ghcr.io/coble-tools/coble:zheng-gen

# Alternatively to run within a docker (eg for a mac) these commands iopen them as containers 
# and mount the current working directory to /workspace in the container. 
# The working directory is set to /workspace so that you can run the scripts from the container.

# docker run  --platform linux/amd64 --rm -it -v .:/workspace -w /workspace ghcr.io/coble-tools/coble:scherer-rep
# docker run  --platform linux/amd64 --rm -it -v .:/workspace -w /workspace ghcr.io/coble-tools/coble:zheng-rep




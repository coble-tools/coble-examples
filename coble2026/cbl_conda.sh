#!/bin/bash

echo "This script is used to build the conda environments for the Coble 2026 paper"

echo "Building conda for Panel(a): Recreate -  Zhang et al (2026) ..."
coble build --recipe Repeat/cbl/recipe.cbl --env zhang-rec --rebuild

echo "Building conda for Panel(b): Reproduce - Scherer et al (2025) ..."
coble build --recipe Reproduce/cbl/recipe.cbl --env scherer-rep --rebuild

echo "Building conda for Panel(c): Replicate - Zheng et al (2017) ..."
coble build --recipe Replicate/cbl/recipe.cbl --env zheng-rep --rebuild

echo "Completed coble builds for Coble 2026 paper"
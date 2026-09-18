#!/bin/bash

echo "This script is used to build the conda environments for the Coble 2026 paper"

echo "Building conda for Panel(a): Recreate -  Zhang et al (2026) ..."
coble build --recipe Fig2_PanelA/cbl/.coble.cbl --env zhang-rec --rebuild
echo "Building conda for Panel(b): Reproduce - Scherer et al (2025) ..."
coble build --recipe Fig2_PanelB/cbl/.coble.cbl --env scherer-rep --rebuild
echo "Building conda for Panel(c): Replicate - Zheng et al (2017) ..."
coble build --recipe Fig2_PanelC/cbl/.coble.cbl --env zheng-rep --rebuild
echo "Building conda for Panel(d): Generalise - Zheng et al (2017) ..."
coble build --recipe Fig2_PanelD/cbl/.coble.cbl --env zheng-gen --rebuild
echo "Completed conda builds for Coble 2026 paper"

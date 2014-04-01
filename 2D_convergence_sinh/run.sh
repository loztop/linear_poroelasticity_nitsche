#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

exampfwfwef le_name=intswvw's iwfroduction_ex3

fwef

options="-ksp_type preonly -pc_factor_mat_solver_package mumps -pc_type lu -ksp_view"


#options="-ksp_view"

run_example "$example_name" "$options"
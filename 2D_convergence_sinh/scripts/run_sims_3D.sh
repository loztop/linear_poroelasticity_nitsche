#!/bin/bash
# My first script

f_prefix="3D_weak_TET"


NE=(3 4 8 10 12) #=NT
DELTA=(0.001)
DELTA_BC=(1000 1000000)

#clpc directory
exe_directory="/users/lorenzb/Dphil/libmesh_git/weak_poro_paper/2D_convergence_sinh/"

#loztop directory
#exe_directory="/home/loztop/Dropbox/Dphil/libmesh_git/weak_poro_paper/2D_convergence_sinh/"

matfiles_dir="data/matfiles/"
data_dir="data/"

res_directory_mat="$exe_directory$matfiles_dir"

res_directory_data="$exe_directory$data_dir"


exe_filename="ex11-opt"


str_nt='NT_'
str_ne='NE_'
str_delta='DELTA_'
str_delta_bc='DELTA_BC_'

for l in ${DELTA_BC[@]}
do
for k in ${DELTA[@]}
do
for i in ${NE[@]}
do





output_file_name_mat="$res_directory_mat$f_prefix"_"$k"_"$l"_"$str_nt$i"_"$str_ne$i"_.mat" "

#output_file_name_data="$res_directory_data$f_prefix"_"$k"_"$str_nt$j"_"$str_ne$i"_" "

#exe_str=" $exe_directory$exe_filename $j $i $output_file_name_mat $output_file_name_data $k" 
exe_str=" $exe_directory$exe_filename $i $i $output_file_name_mat "exudos_dump" $k $l" 

   echo $exe_str

`$exe_str`

done
done
done


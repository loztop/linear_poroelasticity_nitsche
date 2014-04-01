#!/bin/bash
# My first script

f_prefix="conv_linear3DTime_p1p1p0_1stab"

NE=(10 12 18 20 22 24 26)
NT=(10 12 18 20 22 24 26)

#top directory
exe_directory="~/Dropbox/Dphil/libmesh_git/weak_poro_paper/2D_convergence/"

res_directory="~/Dropbox/Dphil/libmesh_git/weak_poro_paper/2D_convergence/data/matfiles/"


exe_filename="ex11-opt"


str_nt='NT_'
str_ne='NE_'

for i in ${NE[@]}
do
for j in ${NT[@]}
do

output_file_name="$res_directory$f_prefix"_"$str_nt$j"_"$str_ne$i"_.mat" "

exe_str="$exe_directory$exe_filename $j $i $output_file_name"
   echo $exe_str

`$exe_str`

done
done




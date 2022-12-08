#!/bin/bash

l=$1
p=$2
r=$3
m=$4


feature_dir="../features/noisy_downsample"
proj_dir="../proj_images/noisy_downsample"
option_dir="./options_downsample_noisy"
init_dir="../mom_lsq_init/noisy_downsample"
results_dir="./results_noisy_downsample"


# create the customized option file
# modify the corresponding base option file options_ptn_"$l"_base accordingly
cat "$option_dir"/options_ptn_"$l"_base > "$option_dir"/options_ptn_"$l"_"$p"_"$m"
#echo "radial_weight $p" >> "$option_dir"/options_ptn_"$l"_"$p"_"$m"
echo "proj_weight $p" >> "$option_dir"/options_ptn_"$l"_"$p"_"$m"

num_proj=`cat $proj_dir/proj_image_2d_out_$m | wc -l`
echo "num_proj $num_proj" >> "$option_dir"/options_ptn_"$l"_"$p"_"$m"

# perform reconstruction
./main "$feature_dir"/spatial_radial_lgwt_feat_noisy "$feature_dir"/spatial_radial_lgwt_map "$feature_dir"/spatial_radial_lgwt "$proj_dir"/proj_image_2d_out_"$m" "$feature_dir"/fourier_ac_decomposition_feat_noisy_"$l" "$feature_dir"/fourier_radial_uniform_map "$feature_dir"/l_idx_fourier_ac_decomposition_"$l" "$option_dir"/options_ptn_"$l"_"$p"_"$m" "$init_dir"/mom_lsq_init_den_3d_"$m" "$results_dir"/output_"$l"_"$p"_"$r"_"$m" "$results_dir"/obj_"$l"_"$p"_"$r"_"$m"



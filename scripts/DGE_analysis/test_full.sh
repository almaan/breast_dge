count_dir="/home/alma/ST-2018/CNNp/data/pre-tile-extraction"
feature_dir="/home/alma/ST-2018/CNNp/DGE/data/YAO_feature_files"
select_for="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/select.samples.dat"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/design.test2.dat"
#output_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res"

for ii in {1..5}; do
    output_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res/"
    Rscript ./DGE_analysis.r --count_dir $count_dir --feature_dir $feature_dir --design_file $design_file --feature_name tumor --select_for $select_for --output_dir ${output_dir}"DGE_res"${ii}".csv" -a --k_members 250 --max_dist -1 --n_samples 4;
done;

count_dir="/home/alma/ST-2018/CNNp/data/pre-tile-extraction"
feature_dir="/home/alma/ST-2018/CNNp/DGE/data/YAO_feature_files"
select_for="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/select.samples2.dat"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/design.test2.dat"
#output_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res"

for ii in {1..3}; do
    output_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res/"
    Rscript ./DGE_analysis.r --count_dir $count_dir --feature_dir $feature_dir --design_file $design_file --feature_name tumor --select_for $select_for --output_dir ${output_dir}"DGE_res"${ii}".csv" -a --k_members 150 --max_dist 3 --n_samples 10;
cat ${output_dir}"DGE_res"${ii}".csv" | sed 's/,/\t/g' | sort -k7 | cut -f1 | head -n 500 | sort -k1 > ${output_dir}"top500_res"${ii}".txt"
done;

count_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression"
feature_dir="/home/alma/ST-2018/CNNp/DGE/data/YAO_feature_files"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_data/design_file_test.dat"


Rscript ./dge_se.r --count_dir $count_dir --feature_dir $feature_dir --design_file $design_file

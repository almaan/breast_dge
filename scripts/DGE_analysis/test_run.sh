count_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/cnt"
feature_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/ft"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_data/design_file_test.dat"


Rscript ./DGE_analysis.r --count_dir $count_dir --feature_dir $feature_dir --design_file $design_file --feature_name tumor

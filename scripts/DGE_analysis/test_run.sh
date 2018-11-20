feature_dir="/home/alma/ST-2018/CNNp/DGE/data/test_data"
count_dir="/home/alma/ST-2018/CNNp/DGE/data/test_data"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_data/design_file_test.dat"
gene_file="/home/alma/ST-2018/CNNp/DGE/data/test_data/PAM50.ENSG.dat"

Rscript ./dge_trial_test.r --count_dir $count_dir --feature_dir $feature_dir --design $design_file --gene_file $gene_file --workers 2

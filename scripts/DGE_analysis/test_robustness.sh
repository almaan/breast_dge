count_dir="/home/alma/ST-2018/CNNp/data/pre-tile-extraction"
feature_dir="/home/alma/ST-2018/CNNp/DGE/data/YAO_feature_files"
select_for="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/robustness"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/robustness/design.in_patient.dat"
patient_list="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/robustness/for.robustness.dat"
output_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res/robustness4"


while read -r patient_id; do
    p_odir=${output_dir}/$patient_id
    if [ ! -d $p_odir ]; then
        mkdir $p_odir;
    fi;
    echo $patient_id > $select_for/"tmp.dat"
    for n in {1..3}; do
        Rscript ./DGE_analysis.r --count_dir $count_dir --feature_dir $feature_dir --design_file $design_file --feature_name tumor --select_for $select_for"/tmp.dat" --output_dir ${p_odir}"/DGE_res"${n}"_"${patient_id}".csv" -a --k_members 100 --max_dist -1 --n_samples 20;
    done;

done<$patient_list
rm $select_for/"tmp.dat"

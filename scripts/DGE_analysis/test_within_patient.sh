count_dir="/home/alma/ST-2018/CNNp/data/pre-tile-extraction"
feature_dir="/home/alma/ST-2018/CNNp/DGE/data/YAO_feature_files"
select_for="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/in_patient"
design_file="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/in_patient/design.in_patient.dat"
patient_list="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/env/in_patient/patient.with.two.dat"
output_dir="/home/alma/ST-2018/CNNp/DGE/data/test_dge_expression/res/in_patient"


while read -r patient_id; do
    p_odir=${output_dir}/$patient_id
    if [ ! -d $p_odir ]; then
        mkdir $p_odir;
    fi;
    echo $patient_id > $select_for/"tmp.dat"

    Rscript ./DGE_analysis.r --count_dir $count_dir --feature_dir $feature_dir --design_file $design_file --feature_name tumor --select_for $select_for"/tmp.dat" --output_dir ${p_odir}"/DGE_res"${patient_id}".csv" -a --k_members 200 --max_dist 3 --n_samples 50;

done<$patient_list
rm $select_for/"tmp.dat"

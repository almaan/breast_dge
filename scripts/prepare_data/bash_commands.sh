#move .gz files to normal .tsv files

for ii in *.gz; do name=$( echo $ii |  egrep "^(.*?)\.tsv" -o); mv $ii $name; done;



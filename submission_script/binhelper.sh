#!/usr/bin/bash

manifest_file=$1
out_file=$(basename $manifest_file).log
sample_id=$(basename $manifest_file |  rev | cut -f2- -d"_" | rev)
ass=`echo $sample_id | cut -f1 -d_`
fold=bins/$ass/${sample_id}

tar xf /home/moritz/data/data_submit/${fold}.tar.gz $fold/${sample_id}.fna.gz

java -jar ~/share/webin-cli-3.0.0.jar -context=genome -manifest=$manifest_file -userName='Webin-41995' -password='775632' -submit > $out_file

run_id=$(grep "The submission has been completed successfully. The following analysis" $out_file | rev | cut -f1 -d" "| rev)

if grep -q "The submission has been completed successfully." $out_file ;
then
  echo $sample_id $run_id >> bin_completed_submissions.tsv
  mv $manifest_file data/bin_manifests_done
  mv $out_file data/bin_manifests_done
  echo $1 done successfully
else
  mv $out_file data/bin_manifests_fail
#  echo $1 failed
fi

rm $fold/${sample_id}.fna.gz

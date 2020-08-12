#!/usr/bin/bash

manifest_file=$1
out_file=$(basename $manifest_file).log
sample_id=$(basename $manifest_file |  rev | cut -f2- -d"_" | rev)

java -jar ~/share/webin-cli-3.0.0.jar  -context=genome -manifest=$manifest_file -userName="Webin-41995" -password="775632" -submit > $out_file

run_id=$(grep "The submission has been completed successfully. The following analysis" $out_file | rev | cut -f1 -d" "| rev)

if grep -q "The submission has been completed successfully." $out_file ;
then
  echo $sample_id $exp_id $run_id >> sag_completed_submissions.tsv
  mv $manifest_file data/manifests_done
  mv $out_file data/manifests_done
else
  mv $out_file data/manifests_fail
fi

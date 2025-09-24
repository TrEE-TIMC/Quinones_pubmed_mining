wget ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/*

for element in *xml.gz ;
  do md5_file=$(md5sum $element | awk '{print $1}');
  md5_tocompare=$(cat $element.md5 | awk '{print $2}');
  if [ $md5_file != $md5_tocompare ]; then
    echo $element" is uncomplete"
  fi;
done;

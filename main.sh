### 1. ensure conservation of miRNA between rats and humans
# 1. mirdb scrape
#mirna="miR-137-5p"
cat mirnas.list | while read mirna;do
  python mirdb_scrape.py rno-"$mirna"
  python mirdb_scrape.py hsa-"$mirna"
done
# 2. calc containment
echo "cutoff","mirna","rno","hsa","shared_hits","jaccard%","containment%","wt_containment" > containment.report
cat mirnas.list | while read mirna;do
  for cutoff in {95..50..5};do echo $cutoff;
    tail -n+2 rno-"$mirna".tab | awk -F"\t" -v ct=$cutoff '{if($3>=ct)print toupper($5)}' | sed 's/ //g' | sort | uniq > rno_list
    tail -n+2 hsa-"$mirna".tab | awk -F"\t" -v ct=$cutoff '{if($3>=ct)print toupper($5)}' | sed 's/ //g' | sort | uniq > hsa_list
    rno=$(cat rno_list | wc -l)
    hsa=$(cat hsa_list | wc -l)
    if [ $rno -lt $hsa ];then sm=$rno;else sm=$hsa;fi
    total=$((rno+hsa))
    hits=$(comm -12 rno_list hsa_list | wc -l)
    jac=0
    con=0
    wcon=0
    if [ $total -gt 0 ];then jac=$(echo "scale=4; a = $hits / $total; a * 100" | bc | sed 's/00$/%/'); fi
    if [ $sm -gt 0 ];then con=$(echo "scale=4; a = $hits / $sm; a * 100" | bc | sed 's/00$/%/');
                          wcon=$(echo "scale=4; a = $hits / $sm; a * $hits" | bc | sed 's/00$//');fi
    echo "$cutoff","$mirna","$rno","$hsa","$hits","$jac","$con","$wcon" >> containment.report
  done
done

grep ^$cutoff containment.report | cut -d"," -f2 > wcon.csv
grep ^$cutoff containment.report | cut -d"," -f2 > jacc.csv
for cutoff in {95..50..5};do echo $cutoff;
  grep ^$cutoff containment.report | cut -d"," -f8 > temp
  cp wcon.csv wcon.csv_temp
  paste -d"," wcon.csv_temp temp > wcon.csv
  grep ^$cutoff containment.report | cut -d"," -f6 > temp
  cp jacc.csv jacc.csv_temp
  paste -d"," jacc.csv_temp temp > jacc.csv
done
rm temp wcon.csv_temp jacc.csv_temp


#### Assess disease impact
wget --quiet http://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz
gunzip all_gene_disease_associations.tsv.gz



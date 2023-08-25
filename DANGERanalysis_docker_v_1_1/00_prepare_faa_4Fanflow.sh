# prepare protein fasta file for Fanflow
# in *pep.all.fa (Ensembl), there are several entries without "gene_symbol"
# this script make it insert "gene_symbol" content into entries without it

db_file=($1 \
)

log=prepare.fasta.error.log

# cd ./$db_dir &&
for i in ${db_file[@]};
do
mkdir ./prepare.faa_${i}/ &&
echo run_00_prepare_faa_4Fanflow.sh_${i}_log >> $log &&
seqkit grep -n -r -p gene_symbol ./${i} \
-o ./prepare.faa_${i}/wGeneSymbol.fa 2>> $log &&
seqkit grep -n -r -v -p gene_symbol ./${i} \
-o ./prepare.faa_${i}/woGeneSymbol.fa 2>> $log &&
seqkit grep -n -r -v -p description ./prepare.faa_${i}/woGeneSymbol.fa \
-o ./prepare.faa_${i}/woGeneSymbol_woDescription.fa 2>> $log &&
seqkit grep -n -r -p description ./prepare.faa_${i}/woGeneSymbol.fa \
-o ./prepare.faa_${i}/woGeneSymbol_wDescription.fa 2>> $log &&
seqkit grep -n -r -v -p description ./prepare.faa_${i}/wGeneSymbol.fa \
-o ./prepare.faa_${i}/wGeneSymbol_woDescription.fa 2>> $log &&
seqkit grep -n -r -p description: ./prepare.faa_${i}/wGeneSymbol.fa \
-o ./prepare.faa_${i}/wGeneSymbol_wDescription.fa 2>> $log &&
seqkit replace -p " description:" -r " gene_symbol:non-available description:" ./prepare.faa_${i}/woGeneSymbol_wDescription.fa \
-o ./prepare.faa_${i}/addGeneSymbol_wDescription.fa 2>> $log &&
seqkit replace -p "\$" -r " description:non-available" ./prepare.faa_${i}/wGeneSymbol_woDescription.fa \
-o ./prepare.faa_${i}/wGeneSymbol_addDescription.fa 2>> $log &&
seqkit replace -p "\$" -r " gene_symbol:non-available description:non-available" ./prepare.faa_${i}/woGeneSymbol_woDescription.fa \
-o ./prepare.faa_${i}/addGeneSymbol_addDescription.fa 2>> $log &&
cd ./prepare.faa_${i}/ &&
cat wGeneSymbol_wDescription.fa wGeneSymbol_addDescription.fa addGeneSymbol_wDescription.fa  addGeneSymbol_addDescription.fa > ../${i}2 2>> ../$log #&&
# cd ../ 
done
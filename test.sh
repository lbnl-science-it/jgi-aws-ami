#! /bin/bash
source .elb-venv/bin/activate
cd nvme_ssd
wget -O TOBG-GENOMES.tar.gz https://figshare.com/ndownloader/files/8849371
tar xf TOBG-GENOMES.tar.gz TOBG_NP-110.fna
../tools/mastiff -o matches.csv  TOBG_NP-110.fna
cat matches.csv | awk -F , '{ if ($2 > "0.99") print $1}'| sed 1d > extract.csv
head extract.csv -n 2 > extract_less.csv
mkdir -p data && cd data
cp ../TOBG_NP-110.fna ./
makeblastdb -in TOBG_NP-110.fna -input_type fasta -dbtype nucl -parse_seqids
aws s3 cp TOBG_NP-110.fna-nucl-metadata.json $S3URL/custom_blastdb/
for f in TOBG_NP-110.fna.n* ; do aws s3 cp $f $S3URL/custom_blastdb/; done
mkdir sra-cache
cat ../extract_less.csv | xargs -I{} sh -c ' time aws s3 cp s3://sra-pub-run-odp/sra/{}/{} --no-sign-request ./sra-cache '
cat ../extract_less.csv | xargs -I{} sh -c ' time fasterq-dump -e 32 --fasta-unsorted --skip-technical -t ./sra-cache ./sra-cache/{} -Z | aws s3 cp - $S3URL/queries/{}.fa --only-show-errors '
parallel -t --jobs 2 --tag -q -a ../extract_less.csv elastic-blast submit --dry-run --query $S3URL/queries/{}.fa --db $S3URL/custom_blastdb/TOBG_NP-110.fna --program blastn --num-nodes 8 --batch-len 1000000000 --results $S3URL/results/results_parallel/output_{} -- -task megablast -word_size 28 -evalue 0.00001 -max_target_seqs 10000000 -perc_identity 97 -outfmt "6 std qlen slen qcovs"  -mt_mode 1

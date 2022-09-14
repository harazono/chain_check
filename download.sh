wget -c -nc -O hg19_last_modified_2020-01-17.fa.gz               https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/latest/hg19.fa.gz
wget -c -nc -O hg38_last_modified_2020-03-16.fa.gz               https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
wget -c -nc -O hg19ToHg38_last_modified_2013-12-31.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19_last_modified_2020-01-17.fa.gz hg38_last_modified_2020-03-16.fa.gz hg19ToHg38_last_modified_2013-12-31.over.chain.gz
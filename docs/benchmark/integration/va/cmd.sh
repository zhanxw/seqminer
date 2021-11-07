# cp -sf ../../data/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz* .
xtime ~/bin/Rscript annotate.R ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz anno.chr1.chunk1000.vcf 1000 2>&1 |tee annotate.chr1.chunk.1000.timing &
xtime ~/bin/Rscript annotate.R ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz anno.chr1.chunk5000.vcf 5000 2>&1 |tee annotate.chr1.chunk.5000.timing &
xtime ~/bin/Rscript annotate.R ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz anno.chr1.chunk10000.vcf 10000 2>&1 |tee annotate.chr1.chunk.10000.timing &
xtime ~/bin/Rscript annotate.R ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz anno.chr1.chunk20000.vcf 20000 2>&1 |tee annotate.chr1.chunk.20000.timing &

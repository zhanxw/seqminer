java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.66 ../data/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
bgzip snpEff.annotated.wgs.vcf

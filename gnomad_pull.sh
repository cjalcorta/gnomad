#!/bin/bash

SECONDS=0

# When run, the user must enter either 'IDUA' or 'GAA' as a command line argument. This ultimately informs what data is downloaded
if [[ $1 = 'IDUA' ]]
then
	genome_url='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.4.vcf.bgz 4:980285-998817'
	exome_url='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.4.vcf.bgz 4:980285-998817'
elif [[ $1 = 'GAA' ]]
then
	genome_url='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.17.vcf.bgz 17:78074855-78094179'
	exome_url='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.17.vcf.bgz 17:78074855-78094179'
else
	echo "ERROR: Please choose either IDUA or GAA"
	exit 1
fi

# tabix is used to download the VCF files from gnomad through google
echo "Downloading Genome Data..."
tabix -h $genome_url,  > genome.vcf

echo "Downloading Exome Data..."
tabix -h $exome_url > exome.vcf

# the VCF files are compressed so that they will be compatibile with bcftools
echo "Compressing Genome Data..."
bgzip genome.vcf

echo "Compressing Exome Data..."
bgzip exome.vcf

# Because the bcftools merge function drops AC and AN by default (and there is no simple way to override this), these columns are renamed to 'allelecount' and 'allelenumber' respectively.
echo "Re-tagging Exome Data..."
bcftools view -O v exome.vcf.gz \
  | sed -e 's/\([;=[:space:]]\)AC\([,;=[:space:]]\)/\1allelecount\2/' \
        -e 's/\([;=[:space:]]\)AN\([,;=[:space:]]\)/\1allelenumber\2/' \
  | bcftools view -O z -o exome.re-tagged.vcf.gz 2>/dev/null

echo "Re-tagging Genome Data..."
bcftools view -O v genome.vcf.gz \
  | sed -e 's/\([;=[:space:]]\)AC\([,;=[:space:]]\)/\1allelecount\2/' \
        -e 's/\([;=[:space:]]\)AN\([,;=[:space:]]\)/\1allelenumber\2/' \
  | bcftools view -O z -o genome.re-tagged.vcf.gz 2>/dev/null

# Several files no longer needed are removed
rm -f genome.vcf.gz exome.vcf.gz

if [[ $1 = "IDUA" ]]
then
	rm -f gnomad.genomes.r2.1.1.sites.4.vcf.bgz.tbi gnomad.exomes.r2.1.1.sites.4.vcf.bgz.tbi
else
	rm -f gnomad.genomes.r2.1.1.sites.17.vcf.bgz.tbi gnomad.exomes.r2.1.1.sites.17.vcf.bgz.tbi
fi

# Since the data has been filtered for specific ranges, it must now be re-indexed with tabix before merge
echo "Indexing Genome Data..."
tabix -p vcf genome.re-tagged.vcf.gz
echo "Indexing Exome Data..."
tabix -p vcf exome.re-tagged.vcf.gz

# The data is merged using bcftools merge function. The --info-rules parameter tells the function to sum all of the listed fields (by default it would choose the first one)
# The -m none parameter makes it so that multiallelic variants are treated as separate
echo "Merging Data..."
bcftools merge exome.re-tagged.vcf.gz genome.re-tagged.vcf.gz -m none --info-rules 'allelecount:sum','allelenumber:sum','AC_afr:sum','AN_afr:sum','nhomalt_afr:sum','AC_amr:sum','AN_amr:sum','nhomalt_amr:sum','AC_eas:sum','AN_eas:sum','nhomalt_eas:sum','AC_nfe:sum','AN_nfe:sum','nhomalt_nfe:sum','AC_fin:sum','AN_fin:sum','nhomalt_fin:sum','AC_asj:sum','AN_asj:sum','nhomalt_asj:sum','AC_sas:sum','AN_sas:sum','nhomalt_sas:sum','AC_oth:sum','AN_oth:sum','nhomalt_oth:sum' -o full.vcf.gz 2>/dev/null
rm -f exome.re-tagged.vcf.gz genome.re-tagged.vcf.gz exome.re-tagged.vcf.gz.tbi genome.re-tagged.vcf.gz.tbi

# Now selecting only the fields that we are interested in, and converting the file to CSV.
echo "Filtering and converting data to CSV..."
bcftools query --print-header -f '%CHROM,%POS,%REF,%ALT,%allelecount,%allelenumber,%AC_afr,%AN_afr,%nhomalt_afr,%AC_amr,%AN_amr,%nhomalt_amr,%AC_eas,%AN_eas,%nhomalt_eas,%AC_nfe,%AN_nfe,%nhomalt_nfe,%AC_fin,%AN_fin,%nhomalt_fin,%AC_asj,%AN_asj,%nhomalt_asj,%AC_sas,%AN_sas,%nhomalt_sas,%AC_oth,%AN_oth,%nhomalt_oth\n' 'full.vcf.gz' | tr "\t" "," > gnomad_pull_$1.csv
rm -f full.vcf.gz

#awk -F , -v OFS=, '$7=$5/$6' out.csv > finalout.csv
#echo "Calculated Total Allele Frequency"INFO/
#rm -f out.csv

duration=$SECONDS
echo "Done, $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
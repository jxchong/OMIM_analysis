set -eo pipefail


DATESTRING=$1;


perl code/parse_morbidmap.pl --in raw_download_${DATESTRING}/morbidmap.txt --out ${DATESTRING}.morbidmap.parsed.txt
perl code/parse_omimtxtZ_count_NGS_year.pl --in raw_download_${DATESTRING}/omim.txt.gz --out ${DATESTRING}.omimtxtZ.parsed.NGS.year.txt
perl code/combine_omimtxtZ_morbidmap_count_NGS_year.pl --morbidmap ${DATESTRING}.morbidmap.parsed.txt --mim2gene raw_download_${DATESTRING}/mim2gene.txt --omimtxt ${DATESTRING}.omimtxtZ.parsed.NGS.year.txt --out ${DATESTRING}.combinedOMIM.mentionsNGS.year.txt

Rscript --vanilla code/analyze_omimtxtZ_morbidmap.R ${DATESTRING}
Rscript --vanilla code/NGS_gene_discoveries_by_year.R ${DATESTRING}


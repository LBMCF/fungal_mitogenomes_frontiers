#This script is part of supplementary documents of "Global Characterization of Fungal Mitogenomes: New Insights on Genomic Diversity and Dynamism of Coding Genes and Accessory Elements"
#submitted to Frontiers in Microbiology, research topic "Mitochondrial Genomes and Mitochondrion Related Gene Insights to Fungal Evolution"
#DOI: 10.3389/fmicb.2021.787283
#Authors:  Paula L. C. Fonseca, Ruth B. De-Paula, Daniel S. Araújo, Luiz Marcelo R. Tomé, Thairine Mendes-Pereira, Wenderson F. Rodrigues, 
#Luiz-Eduardo Del-Bem, Eric R. Aguiar, and Aristóteles Góes-Neto

#Bash script that uses mfannot files(.tbl) to generate 8 output files: one containing the number of features per category (uorfs, introns, HEs/HEGs, core genes, 
#trna, dpo+rpo+rps3) (_counts_final); six containing the coordinates and sizes for each feature in each category (_coords_size_ok); and one that contains all _coords_size_ok files together
#(_all_features_size), to facilitate running feature overlap python script.
#
#The files used are as follow: 
# Insert all .tbl annotation files in the working folder, containing the .tbl extension.

#******************************************************************************#
#                             Run the code in Bash                             #
#******************************************************************************#

for j in *.tbl ; do

# uorfs
    grep -i orf ${j} | grep gene | cut -f 5 > list_orf
    grep -i -B 1 lagli ${j} | grep -i orf | sed 's/^.*G-//g' > list_he_orf
    grep -i -B 1 giy ${j} | grep -i orf | sed 's/^.*G-//g' >> list_he_orf
    grep -i -B 1 hnh ${j} | grep -i orf | sed 's/^.*G-//g' >> list_he_orf
    grep -w -v -f list_he_orf list_orf | sed -e 's/\r//g' > list_uorf
    
    for i in $(cat list_uorf) ; do
        grep -i -B 2 $i ${j} | grep gene | cut -f 1,2 | gawk 'NF' > coord && sed "s/^/${i}\t/g" coord >> uorf_coords && rm coord ; 
    done
    
# introns
    grep -i -A 2 intron ${j} | cut -f 1,2 | gawk 'NF' > intron_coord
    grep -i -A 2 intron ${j} | cut -f 5 | grep Group | sed 's/\r/\n/g' | gawk 'NF' > intron_name
    paste intron_name intron_coord | sed 's/ /_/g' > intron_coords
    rm intron_coord && rm intron_name
    
# HEGs
    grep -i -B 3 lagli ${j} | grep CDS | cut -f 1,2 | gawk 'NF' > coord && sed "s/^/LAGLIDADG\t/g" coord > he_coords && rm coord
    grep -i -B 3 giy ${j} | grep CDS | cut -f 1,2 | gawk 'NF' > coord && sed "s/^/GIY\t/g" coord >> he_coords && rm coord
    grep -i -B 3 hnh ${j} | grep CDS | cut -f 1,2 | gawk 'NF' > coord && sed "s/^/HNH\t/g" coord >> he_coords && rm coord
    
# tRNAs
    grep -i trna ${j} | grep -P "\ttRNA" | cut -f 1,2 > trna_coord
    grep -i -A 2 trna ${j} | grep "G-" | sed 's/^.*G-//g' | sed -e 's/\r//g' > trna_name
    paste trna_name trna_coord > trna_coords
    rm trna_coord && rm trna_name

# core genes
    grep -i gene ${j} | grep -i -v trn | grep -i -v dpo | grep -i -v rpo  | grep -i -v rps3 | grep -i -v orf | grep -i -v rn | cut -f 5 | sed '/^\s*$/d' | sed -e 's/\r//g' > core_gene_list
    
    for i in $(cat core_gene_list) ; do
        grep -i -B 1 -P "${i}" ${j} | cut -f 1,2 | sed '/^\s*$/d' | grep -v "-" >> core_gene_coord ;
    done

    paste core_gene_list core_gene_coord > core_gene_coords
    rm core_gene_list && rm core_gene_coord

# dpo, rpo and rps3 genes
    grep -i -B 1 -P "\tdpo" ${j} | cut -f 1,2 | sed '/^\s*$/d' | sed "s/^/dpo\t/g" > dpo_rpo_rps3_coords
    grep -i -B 1 -P "\trpo" ${j} | cut -f 1,2 | sed '/^\s*$/d' | sed "s/^/rpo\t/g" >> dpo_rpo_rps3_coords
    grep -i -B 1 -P "\trps3" ${j} | cut -f 1,2 | sed '/^\s*$/d' | sed "s/^/rps3\t/g" >> dpo_rpo_rps3_coords

# final files
    echo uorf intron he trna core_gene dpo_rpo_rps3 > list_lists

    for i in $(cat list_lists) ; do
       gawk '{print $1 "\t" $2 "\t" $3 "\t" $3-$2}' ${i}_coords > ${i}_coords_size
       grep -P "\t-" ${i}_coords_size | gawk '{print $1 ":" $2 "_" $3 "\t" (-$4+1) "\t" "-"}' > ${i}_coords_size_-
       grep -v -P "\t-" ${i}_coords_size | gawk '{print $1 ":" $2 "_" $3 "\t" ($4+1) "\t" "+"}' > ${i}_coords_size_+
       cat ${i}_coords_size_+ ${i}_coords_size_- > ${j}_${i}_coords_size_ok
       rm ${i}_coords_size && rm ${i}_coords_size_+ && rm ${i}_coords_size_-
       wc -l ${j}_${i}_coords_size_ok > ${j}_${i}_counts ;
    done

    paste ${j}_uorf_counts ${j}_intron_counts ${j}_he_counts ${j}_trna_counts ${j}_core_gene_counts ${j}_dpo_rpo_rps3_counts > ${j}_counts_final
    cat ${j}_uorf_coords_size_ok ${j}_intron_coords_size_ok ${j}_he_coords_size_ok ${j}_trna_coords_size_ok ${j}_core_gene_coords_size_ok ${j}_dpo_rpo_rps3_coords_size_ok > ${j}_all_features_size
    
    rm *list*
    rm *_coords

    # user message
    echo "- Number of uORFs:" &
    more ${j}_uorf_counts
    echo "."
    echo "."
    echo "."
    echo "- Number of introns:" &
    more ${j}_intron_counts
    echo "."
    echo "."
    echo "."
    echo "- Number of HEs:" &
    more ${j}_uorf_counts
    echo "."
    echo "."
    echo "."
    echo "- Number of tRNAs:" &
    more ${j}_trna_counts
    echo "."
    echo "."
    echo "."
    echo "- Number of core genes:" &
    more ${j}_core_gene_counts
    echo "."
    echo "."
    echo "."
    echo "- Number of dpo, rpo and rps3:" &
    more ${j}_dpo_rpo_rps3_counts
    echo "."
    echo "."
    echo "."
    echo "- Features size: see _coords_size_ok files"
    echo "..."

    rm *_counts ;

done

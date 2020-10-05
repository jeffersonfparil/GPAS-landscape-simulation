#!/bin/bash
prefix=$1
# prefix=POP_01
touch pop_prefixes.temp.${prefix}
touch pop_size.temp.${prefix}
touch pop_pheno.temp.${prefix}
# random sampling of individuals per population
while read l
do
  # l=$(head -n1 ${prefix}.idx)
  pop=$(echo $l | cut -d' ' -f1)
  echo ONEPOOL_${pop} >> pop_prefixes.temp.${prefix}
  cut -d, -f1 ONEPOOL_${pop}.csv >> pop_size.temp.${prefix}
  cut -d, -f2 ONEPOOL_${pop}.csv >> pop_pheno.temp.${prefix}
done < ${prefix}.idx
paste pop_prefixes.temp.${prefix} pop_size.temp.${prefix} pop_pheno.temp.${prefix} > pop_pooling_info.temp.${prefix}
# build phenotype files and sort pop_pooling_info.temp.${prefix}
Rscript build_pheno_and_sort.r pop_pooling_info.temp.${prefix} ${prefix}
# build sync file
i=$(cut -d' ' -f1 pop_pooling_info.temp.${prefix} | head -n1)
cat ${i}.sync > ${prefix}.sync1.temp.${prefix}
for i in $(cut -d' ' -f1 pop_pooling_info.temp.${prefix} | tail -n+2)
do
  # i=$(cut -d' ' -f1 pop_pooling_info.temp.${prefix} | head -n2)
  cut -f4 ${i}.sync > ${prefix}.sync2.temp.${prefix}
  paste ${prefix}.sync1.temp.${prefix} ${prefix}.sync2.temp.${prefix} > ${prefix}.sync3.temp.${prefix}
  mv ${prefix}.sync3.temp.${prefix} ${prefix}.sync1.temp.${prefix}
done
mv ${prefix}.sync1.temp.${prefix} ${prefix}.sync
rm *.temp.${prefix}

#!/bin/bash

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### SAMPLE EXECUTION 1:
PoolGPAS_InvUra () {
  #######################
  ### PoolGPAS_InvUra ###
  #######################
  DIR=/data/Lolium/Quantitative_Genetics/PoolGPAS_InvUra
  cd ${DIR}
  echo -e "_R1" > ${DIR}/read_grouping.txt
  echo -e "_R2" >> ${DIR}/read_grouping.txt
  RGRP=${DIR}/read_grouping.txt
  echo -e "IF" > ${DIR}/pool_grouping.txt
  echo -e "IG" >> ${DIR}/pool_grouping.txt
  echo -e "IS" >> ${DIR}/pool_grouping.txt
  echo -e "IT" >> ${DIR}/pool_grouping.txt
  echo -e "UF" >> ${DIR}/pool_grouping.txt
  echo -e "UG" >> ${DIR}/pool_grouping.txt
  echo -e "US" >> ${DIR}/pool_grouping.txt
  echo -e "UT" >> ${DIR}/pool_grouping.txt
  PGRP=${DIR}/pool_grouping.txt
  MAPQ=20 # minimum mapping quality at 1% error rate
  BASQ=20 # minimum base calling quality at 1% error rate
  MAF=0.001 # dummy MAF - will iterate across 0.0001, 0.001, and 0.01 below
  DEPTH=10 # dummy DEPTH - will iterate across 1, 10, and 100 below
  REFLINK=http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta
  # POPOOLATION_DIR=/data/Lolium/Softwares/popoolation2_1201
  POPOOLATION_DIR=/data/Lolium/Softwares/popoolation2_1201
  # SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
  SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
  echo -e "IF\t74,74,74,74,74" > ${DIR}/pool_sizes_per_grouping.txt
  echo -e "IG\t84,84,85,84,84" >> ${DIR}/pool_sizes_per_grouping.txt
  echo -e "IS\t92,92,95,92,92" >> ${DIR}/pool_sizes_per_grouping.txt
  echo -e "IT\t89,89,91,89,89" >> ${DIR}/pool_sizes_per_grouping.txt
  echo -e "UF\t66,66,67,66,66" >> ${DIR}/pool_sizes_per_grouping.txt
  echo -e "UG\t71,71,71,71,71" >> ${DIR}/pool_sizes_per_grouping.txt
  echo -e "US\t80,80,84,80,80" >> ${DIR}/pool_sizes_per_grouping.txt
  echo -e "UT\t80,80,84,80,80" >> ${DIR}/pool_sizes_per_grouping.txt
  POOLSIZES=${DIR}/pool_sizes_per_grouping.txt
  echo -e "IF1,0.860945945945946\nIF2,0.888243243243243\nIF3,0.901351351351351\nIF4,0.912162162162162\nIF5,0.928378378378378" > IF_pheno.csv
  echo -e "IG1,0.236309523809524\nIG2,0.313452380952381\nIG3,0.363764705882353\nIG4,0.415357142857143\nIG5,0.553809523809524" > IG_pheno.csv
  echo -e "IS1,0.585978260869565\nIS2,0.68554347826087\nIS3,0.737894736842105\nIS4,0.781739130434783\nIS5,0.825760869565217" > IS_pheno.csv
  echo -e "IT1,0.647078651685393\nIT2,0.718089887640449\nIT3,0.75\nIT4,0.77876404494382\nIT5,0.80876404494382" > IT_pheno.csv
  echo -e 'Pheno_name="IF";\nsig=0.0260509560447872;\nMIN=0.77;\nMAX=0.95;\nperc=[0.2,0.4,0.6,0.8];\nq=[0.88,0.9,0.91,0.92];' > IF_pheno.py
  echo -e 'Pheno_name="IG";\nsig=0.117251938971709;\nMIN=0.14;\nMAX=0.82;\nperc=[0.199524940617577,0.399049881235154,0.600950118764846,0.800475059382423];\nq=[0.28,0.34,0.38,0.45];' > IG_pheno.py
  echo -e 'Pheno_name="IS";\nsig=0.0918815055638145;\nMIN=0.26;\nMAX=0.87;\nperc=[0.198704103671706,0.397408207343413,0.602591792656587,0.801295896328294];\nq=[0.66,0.71,0.76,0.8];' > IS_pheno.py
  echo -e 'Pheno_name="IT";\nsig=0.0593055621975472;\nMIN=0.52;\nMAX=0.86;\nperc=[0.19910514541387,0.398210290827741,0.601789709172259,0.80089485458613];\nq=[0.7,0.74,0.76,0.79];' > IT_pheno.py
  echo -e "UF1,0.873787878787879\nUF2,0.897575757575758\nUF3,0.907164179104478\nUF4,0.917272727272727\nUF5,0.932272727272727" > UF_pheno.csv
  echo -e "UG1,0.194366197183099\nUG2,0.26056338028169\nUG3,0.301549295774648\nUG4,0.351971830985916\nUG5,0.435492957746479" > UG_pheno.csv
  echo -e "US1,0.481\nUS2,0.598125\nUS3,0.659166666666667\nUS4,0.71875\nUS5,0.788875" > US_pheno.csv
  echo -e "UT1,0.61275\nUT2,0.6915\nUT3,0.726428571428571\nUT4,0.753625\nUT5,0.785" > UT_pheno.csv
  echo -e 'Pheno_name="UF";\nsig=0.0256519853478083;\nMIN=0.62;\nMAX=0.99;\nperc=[0.199395770392749,0.398791540785498,0.601208459214502,0.800604229607251];\nq=[0.89,0.9,0.91,0.92];' > UF_pheno.py
  echo -e 'Pheno_name="UG";\nsig=0.0884099841443897;\nMIN=0.05;\nMAX=0.73;\nperc=[0.2,0.4,0.6,0.8];\nq=[0.24,0.28,0.324,0.38];' > UG_pheno.py
  echo -e 'Pheno_name="US";\nsig=0.112436271794714;\nMIN=0.16;\nMAX=0.84;\nperc=[0.198019801980198,0.396039603960396,0.603960396039604,0.801980198019802];\nq=[0.56,0.63,0.69,0.75];' > US_pheno.py
  echo -e 'Pheno_name="UT";\nsig=0.0635954161091974;\nMIN=0.47;\nMAX=0.85;\nperc=[0.198019801980198,0.396039603960396,0.603960396039604,0.801980198019802];\nq=[0.67,0.71,0.74,0.77];' > UT_pheno.py
  GENOME=${DIR}/REF/lope_V1.0.fasta
  # DB=/data/BlastDB/NCBI_NT0060
  DB=/data/BlastDB/NCBI_NT0060
  ### Setup directories, align, pileup, and synchronize, and (DO NOT PERFORM SNP filtering, Pool-GPAS and peak analysis yet)
  SETUP=FALSE     # FALSE TOO SINCE THE LINK TO THE REFERENCE GENOME IS DEAD NOW-2020-01-30 (http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta)
  ALIGN=TRUE      # NOTE: SETUP DIRECTORIES MANUALLY!!! SEE PoolGPAS.sh (i.e. mkdir FASTQ, REF, SAM, and VCF)
  PILEUP=TRUE
  SYNC=TRUE
  FILTER=FALSE    # FALSE
  POOLGPAS=FALSE  # FALSE
  PEAKS=FALSE     # FALSE
  POOLGPCV=FALSE  # FALSE
  ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
  ### Iterate across different MAF and DEPTH filtering parameters and
  ### perform SNP filtering, Pool-GPAS and peak analysis
  for MAF in 0.001 0.0001 0.01
  do
    # MAF=0.001
    for DEPTH in 10 50 100 1
    do
      # DEPTH=100
      ### Filter SNPs by MAF and DEPTH
      SETUP=FALSE
      ALIGN=FALSE
      PILEUP=FALSE
      SYNC=FALSE
      FILTER=TRUE     # FILTER!
      POOLGPAS=FALSE
      PEAKS=FALSE
      POOLGPCV=FALSE
      ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
      ### Pool-GPAS and Pool-GPAS cross-validataion
      SETUP=FALSE
      ALIGN=FALSE
      PILEUP=FALSE
      SYNC=FALSE
      FILTER=FALSE
      POOLGPAS=TRUE   # TRUE
      PEAKS=TRUE      # TRUE
      POOLGPCV=FALSE   # FALSE
      ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
      mv GPAS/ GPAS_MAF${MAF}_DEPTH${DEPTH}/
      # save DEPTH- and MAF- filtered sync, allele frequency, and Fst data
      mv ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync GPAS_MAF${MAF}_DEPTH${DEPTH}/
      rm ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_MAF${MAF}_DEPTH${DEPTH}.sync # clean-up
      mv ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv GPAS_MAF${MAF}_DEPTH${DEPTH}/
      mv ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_COVARIATE_FST.csv GPAS_MAF${MAF}_DEPTH${DEPTH}/
      mkdir GPAS/
    done
  done
}

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### SAMPLE EXECUTION 2:
PoolGPAS_60SEAu () {
  #######################
  ### PoolGPAS_60SEAu ###
  ####################### excluding ACC41, ACC42 and ACC58 because they have less than 10 individual plant representatives (due to low seed germination rates)
  DIR=/data/Lolium/Quantitative_Genetics/PoolGPAS_60SEAu
  cd ${DIR}
  echo -e "_R1" > ${DIR}/read_grouping.txt
  echo -e "_R2" >> ${DIR}/read_grouping.txt
  RGRP=${DIR}/read_grouping.txt
  echo -e "ACC" > ${DIR}/pool_grouping.txt
  PGRP=${DIR}/pool_grouping.txt
  MAPQ=20 # minimum mapping quality at 1% error rate
  BASQ=20 # minimum base calling quality at 1% error rate
  MAF=0.001 # dummy MAF - will iterate across 0.0001, 0.001, and 0.01 below
  DEPTH=10 # dummy DEPTH - will iterate across 1, 10, and 100 below
  REFLINK=http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta
  POPOOLATION_DIR=/data/Lolium/Softwares/popoolation2_1201
  SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
  # echo -e "ACC01,22\nACC02,38\nACC03,41\nACC04,10\nACC05,40\nACC06,21\nACC07,31\nACC08,36\nACC09,29\nACC10,35\nACC11,36\nACC12,38\nACC13,33\nACC14,42\nACC15,39\nACC16,42\nACC17,22\nACC18,39\nACC19,42\nACC20,37\nACC21,42\nACC22,39\nACC23,38\nACC24,20\nACC26,33\nACC27,26\nACC28,42\nACC29,38\nACC30,41\nACC31,39\nACC35,41\nACC32,37\nACC33,12\nACC34,29\nACC36,29\nACC37,29\nACC38,38\nACC39,33\nACC40,42\nACC41,1\nACC42,4\nACC43,21\nACC44,19\nACC45,21\nACC46,23\nACC47,41\nACC48,41\nACC49,40\nACC50,37\nACC51,42\nACC52,40\nACC53,38\nACC54,39\nACC55,22\nACC57,42\nACC58,5\nACC59,39\nACC60,40\nACC61,42\nACC62,42" > pool_sizes.csv
  echo -e "ACC01,22\nACC02,38\nACC03,41\nACC04,10\nACC05,40\nACC06,21\nACC07,31\nACC08,36\nACC09,29\nACC10,35\nACC11,36\nACC12,38\nACC13,33\nACC14,42\nACC15,39\nACC16,42\nACC17,22\nACC18,39\nACC19,42\nACC20,37\nACC21,42\nACC22,39\nACC23,38\nACC24,20\nACC26,33\nACC27,26\nACC28,42\nACC29,38\nACC30,41\nACC31,39\nACC35,41\nACC32,37\nACC33,12\nACC34,29\nACC36,29\nACC37,29\nACC38,38\nACC39,33\nACC40,42\nACC43,21\nACC44,19\nACC45,21\nACC46,23\nACC47,41\nACC48,41\nACC49,40\nACC50,37\nACC51,42\nACC52,40\nACC53,38\nACC54,39\nACC55,22\nACC57,42\nACC59,39\nACC60,40\nACC61,42\nACC62,42" > pool_sizes.csv
  echo -e "ACC" > name.temp
  cut -d, -f2 pool_sizes.csv > size.temp
  sed -zi "s/\n/,/g" size.temp
  paste name.temp size.temp > ${DIR}/pool_sizes_per_grouping.txt
  sed -i 's/,$//g' ${DIR}/pool_sizes_per_grouping.txt
  POOLSIZES=${DIR}/pool_sizes_per_grouping.txt
  GENOME=${DIR}/REF/lope_V1.0.fasta
  DB=/data/BlastDB/NCBI_NT0060
  rm name.temp size.temp
  ### prep phenotypes
  # echo -e "ACC01\nACC02\nACC03\nACC04\nACC05\nACC06\nACC07\nACC08\nACC09\nACC10\nACC11\nACC12\nACC13\nACC14\nACC15\nACC16\nACC17\nACC18\nACC19\nACC20\nACC21\nACC22\nACC23\nACC24\nACC26\nACC27\nACC28\nACC29\nACC30\nACC31\nACC32\nACC33\nACC34\nACC35\nACC36\nACC37\nACC38\nACC39\nACC40\nACC41\nACC42\nACC43\nACC44\nACC45\nACC46\nACC47\nACC48\nACC49\nACC50\nACC51\nACC52\nACC53\nACC54\nACC55\nACC57\nACC58\nACC59\nACC60\nACC61\nACC62" > col1.temp
  # ### batch2 pheno only
  # # echo -e "0\n5.26315789473684\n7.31707317073171\n0\n0\n0\n48.3870967741936\n97.2222222222222\n58.6206896551724\n25.7142857142857\n100\n13.1578947368421\n54.5454545454546\n0\n0\n0\n0\n0\n0\n0\n0\n0\n65.7894736842105\n0\n0\n65.3846153846154\n0\n2.63157894736842\n2.4390243902439\n5.12820512820513\n0\n50\n41.3793103448276\n0\n72.4137931034483\n17.2413793103448\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n10\n0\n0\n12.5\n13.1578947368421\n0\n0\n0\n0\n15.3846153846154\n85\n0\n0" > cleth.temp
  # # echo -e "0\n8.786741713571\n38.0952380952381\n50\n0\n0\n40.9188034188034\n2.38095238095238\n95.1219512195122\n5.12820512820513\n14.2857142857143\n23.0769230769231\n7.69230769230769\n10.8108108108108\n0\n0\n13.6363636363636\n97.5\n2.38095238095238\n0\n0\n0\n10\n100\n16.6666666666667\n21.7391304347826\n32.4324324324324\n7.89473684210526\n0\n26.6341256366723\n0\n5.4054054054054\n0\n16.6666666666667\n18.6206896551724\n16.2162162162162\n0\n0\n9.52380952380952\n100\n100\n66.6666666666667\n0\n4.16666666666667\n3.84615384615385\n90\n31.7073170731707\n16.6666666666667\n85.7142857142857\n11.1111111111111\n2.38095238095238\n5.12820512820513\n0\n0\n7.14285714285714\n0\n19.0476190476191\n9.52380952380952\n10.7142857142857\n74.7967479674797" > glyph.temp
  # # echo -e "0\nNA\n0\n36.3636363636364\nNA\n8.10810810810811\n10\n42.1052631578947\n15.7894736842105\n8.10810810810811\n46.1538461538462\n81.0810810810811\n22.8571428571429\nNA\n13.8888888888889\nNA\n0\nNA\n0\n0\n0\n48.7804878048781\nNA\n13.3333333333333\n6.66666666666667\n23.8095238095238\nNA\n9.09090909090909\n7.31707317073171\n12.1951219512195\n0\n12.5\nNA\n25.6410256410256\nNA\n2.56410256410256\n0\n43.3333333333333\n0\nNA\nNA\nNA\nNA\n25\n9.09090909090909\n10.5263157894737\n26.8292682926829\nNA\n28.2051282051282\nNA\n5\n43.5897435897436\n0\n6.9069069069069\n0\nNA\n41.1764705882353\n5.12820512820513\n35\n20.5882352941176" > sulfo.temp
  # # echo -e "NA\nNA\nNA\n48.1458113526864\nNA\n65\n2.38095238095238\nNA\n50\nNA\n24.4158113526864\n33.7558113526864\n45.8458113526864\nNA\n21.4058113526864\nNA\nNA\nNA\nNA\nNA\n0\n15.7894736842105\nNA\nNA\n2.5\nNA\nNA\nNA\nNA\n18.9189189189189\nNA\n33.3333333333333\nNA\n12.5058113526864\nNA\n44.6558113526864\nNA\nNA\nNA\nNA\nNA\nNA\n15.4858113526864\nNA\nNA\nNA\nNA\n36.3158113526864\n29.7297297297297\nNA\nNA\nNA\n31.5789473684211\n85.7142857142857\nNA\nNA\nNA\n91.8918918918919\n47.3684210526316\n58.5365853658537" > terbu.temp
  # ### merged batch1 and batch2 count data and recalculated resistance levels
  # echo -e "0\n3.03030303030303\n13.2075471698113\n0\n0\n0\n48.3870967741936\n97.2222222222222\n58.6206896551724\n25.7142857142857\n100\n13.1578947368421\n54.5454545454545\n0\n0\n0\n0\n0\n0\n0\n0\n0\n55.1020408163265\n6.89655172413793\n0\n65.3846153846154\n0\n2.63157894736842\n2.17391304347826\n4.8780487804878\n0\n50\n61.3636363636364\n0\n72.5\n17.2413793103448\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n10\n0\n0\n12.5\n13.1578947368421\n0\n3.84615384615385\n0\n0\n15.3846153846154\n85\n0\n0" > cleth.temp
  # echo -e "0\n7.01754385964912\n38.9830508474576\n88.3720930232558\n1.38888888888889\n32.4324324324324\n40\n14\n92.5373134328358\n18.5185185185185\n42.8571428571429\n43.1034482758621\n16.3636363636364\n34.375\n0.826446280991736\n0\n11.7647058823529\n80.5970149253731\n1.47058823529412\n0\n0\n0\n11.7647058823529\n93.6170212765958\n15.7894736842105\n21.7391304347826\n50.9803921568627\n7.14285714285714\n0\n26.4516129032258\n0\n6.25\n0\n20\n15.3846153846154\n17.4603174603175\n0\n0\n11.5384615384615\n100\n100\n66.6666666666667\n27.5\n1.69491525423729\n19.6969696969697\n90.2777777777778\n36.25\n50\n84.375\n19.4444444444444\n2.89855072463768\n4\n6.12244897959184\n0\n5\n0\n21.6981132075472\n9.19540229885057\n20\n79.0476190476191" > glyph.temp
  # echo -e "0\nNA\n0\n36.3636363636364\nNA\n8.10810810810811\n10\n42.1052631578947\n22.0338983050847\n8\n46.1538461538462\n81.0810810810811\n22.8571428571429\n52.1739130434783\n13.8888888888889\nNA\n0\n0\n0\n0\n0\n48.780487804878\nNA\n13.3333333333333\n6.66666666666667\n23.8095238095238\nNA\n9.09090909090909\n7.31707317073171\n12.1951219512195\n0\n12.5\nNA\n25.6410256410256\nNA\n2.56410256410256\n0\n39.1304347826087\n0\nNA\nNA\nNA\nNA\n11.9047619047619\n34.6153846153846\n23.4375\n22.3684210526316\nNA\n28.2051282051282\n56.25\n10.7142857142857\n52\n0\n6.25\n0\nNA\n41.1764705882353\n4.8780487804878\n35\n20.5882352941176" > sulfo.temp
  # echo -e "NA\nNA\nNA\n45.1612903225806\nNA\n64.5161290322581\n2.38095238095238\n100\n50\nNA\n21.4285714285714\n30.7692307692308\n42.8571428571429\nNA\n18.4210526315789\nNA\nNA\nNA\nNA\nNA\n0\n16.3636363636364\nNA\nNA\n4.76190476190476\nNA\nNA\n25\nNA\n18.9189189189189\nNA\n33.3333333333333\nNA\n9.52380952380952\nNA\n41.6666666666667\nNA\nNA\nNA\nNA\nNA\nNA\n12.5\nNA\nNA\nNA\nNA\n33.3333333333333\n33.9285714285714\nNA\nNA\nNA\n29.5454545454545\n85.7142857142857\nNA\nNA\n0\n91.8918918918919\n35.1851851851852\n42.3728813559322" > terbu.temp
  # paste -d, col1.temp cleth.temp > cleth.phen
  # paste -d, col1.temp glyph.temp > glyph.phen
  # paste -d, col1.temp sulfo.temp > sulfo.phen
  # paste -d, col1.temp terbu.temp > terbu.phen
  # rm *.temp
  echo -e 'POP,POOL_SIZE,CLETHODIM,GLYPHOSATE,SULFOMETURON,TERBUTHYLAZINE
  ACC01,22,0,0,0,NA
  ACC02,38,3.03030303030303,7.01754385964912,NA,NA
  ACC03,41,13.2075471698113,38.9830508474576,0,NA
  ACC04,10,0,88.3720930232558,36.3636363636364,45.1612903225806
  ACC05,40,0,1.38888888888889,NA,NA
  ACC06,21,0,32.4324324324324,8.10810810810811,64.5161290322581
  ACC07,31,48.3870967741936,40,10,2.38095238095238
  ACC08,36,97.2222222222222,14,42.1052631578947,100
  ACC09,29,58.6206896551724,92.5373134328358,22.0338983050847,50
  ACC10,35,25.7142857142857,18.5185185185185,8,NA
  ACC11,36,100,42.8571428571429,46.1538461538462,21.4285714285714
  ACC12,38,13.1578947368421,43.1034482758621,81.0810810810811,30.7692307692308
  ACC13,33,54.5454545454545,16.3636363636364,22.8571428571429,42.8571428571429
  ACC14,42,0,34.375,52.1739130434783,NA
  ACC15,39,0,0.826446280991736,13.8888888888889,18.4210526315789
  ACC16,42,0,0,NA,NA
  ACC17,22,0,11.7647058823529,0,NA
  ACC18,39,0,80.5970149253731,0,NA
  ACC19,42,0,1.47058823529412,0,NA
  ACC20,37,0,0,0,NA
  ACC21,42,0,0,0,0
  ACC22,39,0,0,48.780487804878,16.3636363636364
  ACC23,38,55.1020408163265,11.7647058823529,NA,NA
  ACC24,20,6.89655172413793,93.6170212765958,13.3333333333333,NA
  ACC26,33,0,15.7894736842105,6.66666666666667,4.76190476190476
  ACC27,26,65.3846153846154,21.7391304347826,23.8095238095238,NA
  ACC28,42,0,50.9803921568627,NA,NA
  ACC29,38,2.63157894736842,7.14285714285714,9.09090909090909,25
  ACC30,41,2.17391304347826,0,7.31707317073171,NA
  ACC31,39,4.8780487804878,26.4516129032258,12.1951219512195,18.9189189189189
  ACC35,41,0,0,0,NA
  ACC32,37,50,6.25,12.5,33.3333333333333
  ACC33,12,61.3636363636364,0,NA,NA
  ACC34,29,0,20,25.6410256410256,9.52380952380952
  ACC36,29,72.5,15.3846153846154,NA,NA
  ACC37,29,17.2413793103448,17.4603174603175,2.56410256410256,41.6666666666667
  ACC38,38,0,0,0,NA
  ACC39,33,0,0,39.1304347826087,NA
  ACC40,42,0,11.5384615384615,0,NA
  ACC43,21,0,66.6666666666667,NA,NA
  ACC44,19,0,27.5,NA,12.5
  ACC45,21,0,1.69491525423729,11.9047619047619,NA
  ACC46,23,0,19.6969696969697,34.6153846153846,NA
  ACC47,41,0,90.2777777777778,23.4375,NA
  ACC48,41,0,36.25,22.3684210526316,NA
  ACC49,40,10,50,NA,33.3333333333333
  ACC50,37,0,84.375,28.2051282051282,33.9285714285714
  ACC51,42,0,19.4444444444444,56.25,NA
  ACC52,40,12.5,2.89855072463768,10.7142857142857,NA
  ACC53,38,13.1578947368421,4,52,NA
  ACC54,39,0,6.12244897959184,0,29.5454545454545
  ACC55,22,3.84615384615385,0,6.25,85.7142857142857
  ACC57,42,0,5,0,NA
  ACC59,39,15.3846153846154,21.6981132075472,41.1764705882353,0
  ACC60,40,85,9.19540229885057,4.8780487804878,91.8918918918919
  ACC61,42,0,20,35,35.1851851851852
  ACC62,42,0,79.0476190476191,20.5882352941176,42.3728813559322' > phenotype_data.csv
  # ACC41,1,0,100,NA,NA
  # ACC42,4,0,100,NA,NA
  # ACC58,5,0,0,NA,NA
  tail -n+2 phenotype_data.csv | cut -d, -f1 > col1.temp
  tail -n+2 phenotype_data.csv | cut -d, -f3 >  cleth.temp
  tail -n+2 phenotype_data.csv | cut -d, -f4 >  glyph.temp
  tail -n+2 phenotype_data.csv | cut -d, -f5 >  sulfo.temp
  tail -n+2 phenotype_data.csv | cut -d, -f6 >  terbu.temp
  paste -d, col1.temp cleth.temp > cleth.phen
  paste -d, col1.temp glyph.temp > glyph.phen
  paste -d, col1.temp sulfo.temp > sulfo.phen
  paste -d, col1.temp terbu.temp > terbu.phen
  rm *.temp
  ### Setup directories, align, pileup, and synchronize, and (DO NOT PERFORM SNP filtering, Pool-GPAS and peak analysis yet)
  SETUP=FALSE     # FALSE TOO SINCE THE LINK TO THE REFERENCE GENOME IS DEAD NOW-2020-01-30 (http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta)
  ALIGN=TRUE      # NOTE: SETUP DIRECTORIES MANUALLY!!! SEE PoolGPAS.sh (i.e. mkdir FASTQ, REF, SAM, and VCF)
  PILEUP=TRUE
  SYNC=TRUE
  FILTER=FALSE    # FALSE
  POOLGPAS=FALSE  # FALSE
  PEAKS=FALSE     # FALSE
  POOLGPCV=FALSE  # FALSE
  ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
  ### Iterate across different MAF and DEPTH filtering parameters and
  ### perform SNP filtering, Pool-GPAS and peak analysis
  for MAF in 0.001 0.01 0.0001
  do
    # MAF=0.001
    for DEPTH in 10 50 100 1
    do
      # DEPTH=100
      # DEPTH=10
      ### Filter SNPs by MAF and DEPTH
      SETUP=FALSE
      ALIGN=FALSE
      PILEUP=FALSE
      SYNC=FALSE
      FILTER=TRUE     # FILTER!
      POOLGPAS=FALSE
      PEAKS=FALSE
      POOLGPCV=FALSE
      ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
      ### Save a back-up copy to be used across the different herbicide resistance traits
      cp  ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync \
          ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync.bk
      cp  ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv \
          ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv.bk
      ### Iterate across phenotypes and perform Pool-GPAS and peak analysis
      for PHENO in cleth glyph sulfo terbu ### CLETH AND GLYPH DATA TOO BIG FOR 47.2 GB RAM
      # for PHENO in sulfo terbu
      do
        # PHENO=cleth
        echo ${PHENO}
        cp ${PHENO}.phen ACC_pheno.csv
        ### extract phenotype and pool sizes data (remove NA and sort)
        echo -e "
        dat = read.csv(\"ACC_pheno.csv\", header=FALSE)
        sizes = read.csv(\"pool_sizes.csv\", header=FALSE)
        ### remove missing data
        dat = dat[!is.na(dat[,2]), ]
        sizes = sizes[!is.na(dat[,2]), ]
        ### sort by increasing phenotype value and writeout idx sorting to sort the sync and allele frequency files
        idx = order(dat[,2])
        dat = dat[idx, ]
        sizes = sizes[idx, ]
        ### build the python phenotype file
        PHENO = \"Herbicide_resistance\"
        SIG = sd(dat[,2])
        MIN = min(dat[,2])
        MAX = max(dat[,2])
        PERC = cumsum(sizes[,2]/sum(sizes[,2]))[1:(nrow(dat)-1)]
        Q = dat[1:(nrow(dat)-1),2]
        ### writeout
        pool_sizes_per_grouping = data.frame(COL1=\"ACC\", COL2=paste(sizes\$V2, collapse=\",\"))
        OUT = matrix(c( paste0('Pheno_name=\"', PHENO,'\";'),
                        paste0(\"sig=\", SIG, \";\"),
                        paste0(\"MIN=\", MIN, \";\"),
                        paste0(\"MAX=\", MAX, \";\"),
                        paste0(\"perc=[\", paste(PERC, collapse=\",\"), \"];\"),
                        paste0(\"q=[\", paste(Q, collapse=\",\"), \"];\")
                        ), ncol=1)
        write.table(idx, file=\"idx.csv\", sep=\",\", row.names=FALSE, col.names=FALSE, quote=FALSE)
        write.table(dat, file=\"ACC_pheno.csv\", sep=\",\", row.names=FALSE, col.names=FALSE, quote=FALSE)
        write.table(pool_sizes_per_grouping, file=\"pool_sizes_per_grouping.txt\", sep=\"\t\", row.names=FALSE, col.names=FALSE, quote=FALSE)
        write.table(OUT, file=\"ACC_pheno.py\", sep=\"\t\", row.names=FALSE, col.names=FALSE, quote=FALSE)
        " > filter_pheno_data.r
        Rscript filter_pheno_data.r
        ### Extract and sort genotype data from the backup sync and allele frequency files -
        ### based on the populations defined in the idx.csv file where NA were removed and sorted by increasing phenotypic value
        ### SYNC FILE SORTING
        cut -f1-3 ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync.bk > col123.temp
        cut -f4- ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync.bk > col42end.temp
        cut -f$(head -n1 idx.csv) col42end.temp > SORTED.temp
        for j in $(tail -n+2 idx.csv)
        do
          cut -f$j col42end.temp > temp.temp
          paste SORTED.temp temp.temp > sorted.temp
          mv sorted.temp SORTED.temp
        done
        paste col123.temp SORTED.temp > ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync
        rm *.temp
        ### ALLELEFREQ FILE SORTING
        cut -d, -f1-3 ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv.bk > col123.temp
        cut -d, -f4- ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv.bk > col42end.temp
        cut -d, -f$(head -n1 idx.csv) col42end.temp > SORTED.temp
        for j in $(tail -n+2 idx.csv)
        do
          cut -d, -f$j col42end.temp > temp.temp
          paste -d, SORTED.temp temp.temp > sorted.temp
          mv sorted.temp SORTED.temp
        done
        paste -d, col123.temp SORTED.temp > ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv
        rm *.temp
        ### POOL SIZES PER GROUPING FILE SORTING
        mv $POOLSIZES ${POOLSIZES}.bk
        cut -f1 ${POOLSIZES}.bk > col1.temp
        cut -f2- ${POOLSIZES}.bk > col2.temp
        cut -d, -f$(head -n1 idx.csv) col2.temp > SORTED.temp
        for j in $(tail -n+2 idx.csv)
        do
          cut -d, -f$j col2.temp > temp.temp
          paste -d, SORTED.temp temp.temp > sorted.temp
          mv sorted.temp SORTED.temp
        done
        paste col1.temp SORTED.temp > $POOLSIZES
        rm *.temp
        ### Pool-GPAS and Pool-GPAS cross-validataion
        SETUP=FALSE
        ALIGN=FALSE
        PILEUP=FALSE
        SYNC=FALSE
        FILTER=FALSE
        POOLGPAS=TRUE   # TRUE
        PEAKS=TRUE      # TRUE
        POOLGPCV=TRUE   # TRUE
        ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
        mv GPAS/ GPAS_${PHENO}_MAF${MAF}_DEPTH${DEPTH}/
        # save DEPTH- and MAF- filtered sync, allele frequency, and Fst data
        mv ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync GPAS_${PHENO}_MAF${MAF}_DEPTH${DEPTH}/
        rm ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_MAF${MAF}_DEPTH${DEPTH}.sync # clean-up
        mv ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv GPAS_${PHENO}_MAF${MAF}_DEPTH${DEPTH}/
        mv ${DIR}/VCF/ACC_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_COVARIATE_FST.csv GPAS_${PHENO}_MAF${MAF}_DEPTH${DEPTH}/
        mkdir GPAS/
      done
    done
  done
}
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### SAMPLE EXECUTION 3:
PoolGPAS_08SEAu () {
  #######################
  ### PoolGPAS_08SEAu ###
  #######################
  DIR=/data/Lolium/Quantitative_Genetics/PoolGPAS_08SEAu
  cd ${DIR}
  RGRP=${DIR}/read_grouping.txt
  PGRP=${DIR}/pool_grouping.txt
  MAPQ=20 # minimum mapping quality at 1% error rate
  BASQ=20 # minimum base calling quality at 1% error rate
  MAF=0.001 # dummy MAF - will iterate across 0.0001, 0.001, and 0.01 below
  DEPTH=10 # dummy DEPTH - will iterate across 1, 10, and 100 below
  REFLINK=http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta
  POPOOLATION_DIR=/data/Lolium/Softwares/popoolation2_1201
  SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
  POOLSIZES=${DIR}/pool_sizes_per_grouping.txt
  GENOME=${DIR}/REF/lope_V1.0.fasta
  DB=/data/BlastDB/NCBI_NT0060
  ### Setup directories, align, pileup, and synchronize, and (DO NOT PERFORM SNP filtering, Pool-GPAS and peak analysis yet)
  SETUP=FALSE     # FALSE TOO SINCE THE LINK TO THE REFERENCE GENOME IS DEAD NOW-2020-01-30 (http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta)
  ALIGN=TRUE      # NOTE: SETUP DIRECTORIES MANUALLY!!! SEE PoolGPAS.sh (i.e. mkdir FASTQ, REF, SAM, and VCF)
  PILEUP=TRUE
  SYNC=TRUE
  FILTER=FALSE    # FALSE
  POOLGPAS=FALSE  # FALSE
  PEAKS=FALSE     # FALSE
  POOLGPCV=FALSE  # FALSE
  ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
  ### Iterate across different MAF and DEPTH filtering parameters and
  ### perform SNP filtering, Pool-GPAS and peak analysis
  for MAF in 0.001 0.01 0.0001
  do
    # MAF=0.001
    for DEPTH in 10 50 100 1
    do
      # DEPTH=100
      ### Filter SNPs by MAF and DEPTH
      SETUP=FALSE
      ALIGN=FALSE
      PILEUP=FALSE
      SYNC=FALSE
      FILTER=TRUE     # FILTER!
      POOLGPAS=FALSE
      PEAKS=FALSE
      POOLGPCV=FALSE
      ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
      ### Pool-GPAS and Pool-GPAS cross-validataion
      SETUP=FALSE
      ALIGN=FALSE
      PILEUP=FALSE
      SYNC=FALSE
      FILTER=FALSE
      POOLGPAS=TRUE   # TRUE
      PEAKS=TRUE      # TRUE
      POOLGPCV=FALSE   # FALSE
      ${SRC_DIR%/src*}/PoolGPAS.sh ${DIR} ${RGRP} ${PGRP} ${MAPQ} ${BASQ} ${REFLINK} ${POPOOLATION_DIR} ${SRC_DIR} ${MAF} ${DEPTH} ${POOLSIZES} ${GENOME} ${DB} ${SETUP} ${ALIGN} ${PILEUP} ${SYNC} ${FILTER} ${POOLGPAS} ${PEAKS} ${POOLGPCV}
      mv GPAS/ GPAS_MAF${MAF}_DEPTH${DEPTH}/
      # save DEPTH- and MAF- filtered sync, allele frequency, and Fst data
      mv ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync GPAS_MAF${MAF}_DEPTH${DEPTH}/
      rm ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_MAF${MAF}_DEPTH${DEPTH}.sync # clean-up
      mv ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv GPAS_MAF${MAF}_DEPTH${DEPTH}/
      mv ${DIR}/VCF/*_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_COVARIATE_FST.csv GPAS_MAF${MAF}_DEPTH${DEPTH}/
      mkdir GPAS/
    done
  done
}


### Execute which one?
if [ $1 == 1 ]
then
  time PoolGPAS_InvUra
fi

if [ $1 == 2 ]
then
  time PoolGPAS_60SEAu
fi

if [ $1 == 3 ]
then
  time PoolGPAS_08SEAu
fi

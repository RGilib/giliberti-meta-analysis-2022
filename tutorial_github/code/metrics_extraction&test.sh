array1=(
abundance
bool0
bool01
bool001
bool0001
bool00001
)
array2=(
species
genus
order
family
)
array3=(
Random_forest
lsvm
svm
enet
lasso
)
array4=(
auc
f1
precision
recall
auprc
)
array5=(
JieZ_2017
ChngKR_2016
YeZ_2018
RaymondF_2016
QinN_2014
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
LiJ_2017
IjazUZ_2017
NielsenHB_2014
GhensiP_2019_m
GhensiP_2019
Castro-NallarE_2015
Heitz-BuschartA_2016
KosticAD_2015
KarlssonFH_2013
QinJ_2012
HMP_2012
)
array6=(
rf
lsvm
svm
enet
lasso
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
for index5 in ${!array5[*]}
do
(grep ${array4[index4]} ../results/classification/shotgun/${array3[index3]}/${array1[index1]}/${array2[index2]}/${array5[index5]}_${array6[index3]}.txt | head -n1)  >> ../metrics/shotgun/${array3[index3]}/${array1[index1]}/${array2[index2]}/${array4[index4]}.txt
done
done
done
done
done

array1=(
abundance
bool0
)
array2=(
Random_forest
lsvm
svm
enet
lasso
)
array3=(
auc
f1
precision
recall
auprc
)
array4=(
edd_singh
crc_baxter
crc_zeller
crc_wang
crc_chen
ibd_gevers
ibd_morgan
ibd_papa
ibd_willing
hiv_noguerajulian
hiv_dinh
hiv_lozupone
asd_son
asd_kang
t1d_alkanani
t1d_mejialeon
nash_chan
nash_zhu
art_scher
mhe_zhang
par_scheperjans
cdi_schubert
non_cdi_schubert
cdi_vincent
cdi_youngster
ob_goodrich
ob_turnbaug
ob_zupancic
ob_ross
ob_zhu
)
array5=(
rf
lsvm
svm
enet
lasso
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
count=`grep ${array3[index3]} ../results/classification/16s/${array2[index2]}/${array1[index1]}/${array4[index4]}_${array5[index2]}.txt|wc -lc
if [ $count -ne 0 ]
then
(grep ${array3[index3]} ../results/classification/16s/${array2[index2]}/${array1[index1]}/${array4[index4]}_${array5[index2]}.txt | head -n1)  >> ../metrics/16s/${array2[index2]}/${array1[index1]}/${array3[index3]}.txt
else
(echo -n -e ${array4[index4]} "\t" ${array3[index3]} "\t" "\t") >> ../metrics/16s/${array2[index2]}/${array1[index1]}/${array3[index3]}.txt
fi
done
done
done
done

array1=(
Random_forest
lsvm
svm
enet
lasso
)

array2=(
abundance
bool0
bool01
bool001
bool0001
bool00001
)

array3=(
species
genus
order
family
)
array4=(
auc
f1
precision
recall
auprc
)
array5=(
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
)
array6=(
rf
lsvm
svm
enet
lasso
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
for index5 in ${!array5[*]}
do
(grep ${array4[index4]} ../results/classification/shotgun/lodo/${array1[index1]}/${array2[index2]}/${array3[index3]}/${array5[index5]}_crc_${array6[index1]}.txt | head -n1) >> ../results/metrics/shotgun/lodo/${array1[index1]}/${array2[index2]}/${array3[index3]}/${array4[index4]}.txt
done
done
done
done
done

array1=(
bool0
bool01
bool001
bool0001
bool00001
)
array2=(
species
genus
order
family
)
array3=(
Random_forest
lsvm
svm
enet
lasso
)
array4=(
JieZ_2017
ChngKR_2016
YeZ_2018
RaymondF_2016
QinN_2014
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
LiJ_2017
IjazUZ_2017
NielsenHB_2014
GhensiP_2019_m
GhensiP_2019
Castro-NallarE_2015
Heitz-BuschartA_2016
KosticAD_2015
KarlssonFH_2013
QinJ_2012
HMP_2012
)
array5=(
rf
lsvm
svm
enet
lasso
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
python compute_test.py ../results/classification/shotgun/${array3[index3]}/abundance/${array2[index2]}/${array4[index4]}_${array5[index3]}_estimations.txt ../results/classification/shotgun/${array3[index3]}/${array1[index1]}/${array2[index2]}/${array4[index4]}_${array5[index3]}_estimations.txt ../results/test/shotgun/${array3[index3]}/${array1[index1]}/${array2[index2]}/${array4[index4]}_test.txt 
done
done
done
done

array1=(
genus
order
family
)

array2=(
JieZ_2017
ChngKR_2016
YeZ_2018
RaymondF_2016
QinN_2014
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
LiJ_2017
IjazUZ_2017
NielsenHB_2014
GhensiP_2019_m
GhensiP_2019
Castro-NallarE_2015
Heitz-BuschartA_2016
KosticAD_2015
KarlssonFH_2013
QinJ_2012
HMP_2012
)


for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
python compute_test.py ../results/classification/shotgun/Random_forest/abundance/species/${array2[index2]}_rf_estimations.txt ../results/classification/Random_forest/abundance/${array1[index1]}/${array2[index2]}_rf_estimations.txt ../results/test/shotgun/Random_forest/${array1[index1]}/${array2[index2]}_test.txt
done
done


array1=(
Random_forest
)
array2=(
edd_singh
crc_baxter
crc_zeller
crc_wang
crc_chen
ibd_gevers
ibd_morgan
ibd_papa
ibd_willing
hiv_noguerajulian
hiv_dinh
hiv_lozupone
asd_son
asd_kang
t1d_alkanani
t1d_mejialeon
nash_chan
nash_zhu
art_scher
mhe_zhang
par_scheperjans
cdi_schubert
non_cdi_schubert
cdi_vincent
cdi_youngster
ob_goodrich
ob_turnbaug
ob_zupancic
ob_ross
ob_zhu
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
python compute_test.py ../results/classification/16s/${array1[index1]}/abundance/${array2[index2]}_rf.txt ../results/classification/16s/${array1[index1]}/bool0/${array2[index2]}_rf.txt ../results/test/16s/${array1[index1]}/bool0/${array2[index2]}_test.txt
done
done

array1=(
Random_forest
lsvm
svm
enet
lasso
)

array2=(
bool0
bool01
bool001
bool0001
bool00001
)

array3=(
species
genus
order
family
)

array4=(
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
)

for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
python compute_test.py ../results/classification/shotgun/lodo/${array1[index1]}/abundance/${array3[index3]}/${array4[index4]}_crc_${array5[index1]}_estimations.txt ../results/classification/shotgun/lodo/${array1[index1]}/${array2[index2]}/${array3[index3]}/${array4[index4]}_crc_${array5[index1]}_estimations.txt ../results/test/shotgun/lodo/${array1[index1]}/${array2[index2]}/${array3[index3]}/${array4[index4]}_test.txt

done
done
done
done


array1=(
bool0
bool01
bool001
bool0001
bool00001
)
array2=(
species
genus
order
family
)
array3=(
Random_forest
lsvm
svm
enet
lasso
)
array4=(
test auc
test f1
test precision
test recall
)
array5=(
JieZ_2017
ChngKR_2016
YeZ_2018
RaymondF_2016
QinN_2014
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
LiJ_2017
IjazUZ_2017
NielsenHB_2014
GhensiP_2019_m
GhensiP_2019
Castro-NallarE_2015
Heitz-BuschartA_2016
KosticAD_2015
KarlssonFH_2013
QinJ_2012
HMP_2012
)
array6=(
test auc
test f1
test precision
test recall
)
array3=(
Random_forest
lsvm
svm
enet
lasso
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
for index5 in ${!array5[*]}
do
(echo -n -e ${array5[index5]} "\t"; grep ${array4[index4]} ../results/test/shotgun/${array3[index3]}/${array1[index1]}/${array2[index2]}/${array5[index5]}_test.txt | head -n1) >> ../metrics/test/shotgun/${array3[index3]}/${array1[index1]}/${array2[index2]}/${array6[index4]}.txt
done
done
done
done
done

array1=(
genus
order
family
)
array2=(
test auc
test f1
test precision
test recall
)
array3=(
JieZ_2017
ChngKR_2016
YeZ_2018
RaymondF_2016
QinN_2014
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
LiJ_2017
IjazUZ_2017
NielsenHB_2014
GhensiP_2019_m
GhensiP_2019
Castro-NallarE_2015
Heitz-BuschartA_2016
KosticAD_2015
KarlssonFH_2013
QinJ_2012
HMP_2012
)

for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
(grep ${array2[index2]} ../results/test/shotgun/Random_forest/${array1[index1]}/${array3[index3]}_test.txt | head -n1) >> ../metrics/test/shotgun/Random_forest/${array1[index1]}/${array2[index2]}.txt
done
done
Done
array1=(
Random_forest
lsvm
svm
enet
lasso
)
array2=(
edd_singh
crc_baxter
crc_zeller
crc_wang
crc_chen
ibd_gevers
ibd_morgan
ibd_papa
ibd_willing
hiv_noguerajulian
hiv_dinh
hiv_lozupone
asd_son
asd_kang
t1d_alkanani
t1d_mejialeon
nash_chan
nash_zhu
art_scher
mhe_zhang
par_scheperjans
cdi_schubert
non_cdi_schubert
cdi_vincent
cdi_youngster
ob_goodrich
ob_Turnbaug
ob_zupancic
ob_ross
ob_zhu
)
array3=(
test auc
test f1
test precision
test recall
)
array4=(
auc
f1
precision
recall
)
for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
(grep ${array2[index2]} ../results/test/16s/${array1[index1]}/bool0/${array2[index2]}_test.txt | head -n1) >> ../metrics/test/16s/${array1[index1]}/bool0/${array4[index3]}.txt
done
done
done

array1=(
Random_forest
lsvm
svm
enet
lasso
)

array2=(
bool0
bool01
bool001
bool0001
bool00001
)

array3=(
species
genus
order
family
)

array4=(
FengQ_2015
GuptaA_2019
HanniganGD_2017
ThomasAM_2018a
ThomasAM_2018b
VogtmannE_2016
WirbelJ_2018
YachidaS_2019
YuJ_2015
ZellerG_2014
)
array5=(
test auc
test f1
test precision
test recall
)
array6=(
auc
f1
precision
recall
)

for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do
for index4 in ${!array4[*]}
do
for index5 in ${!array5[*]}
do

(grep ${array5[index5]} ../results/test/shotgun/lodo/${array1[index1]}/${array2[index2]}/${array3[index3]}/${array4[index4]}_test.txt | head -n1) >> ../metrics/test/shotgun/lodo/${array1[index1]}/${array2[index2]}/${array3[index3]}/${array6[index5]}.txt
done
done
done
done
done




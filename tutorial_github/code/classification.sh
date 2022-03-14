array1=(
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
study_condition
body_site
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
ACVD
AD
BD
cephalosporins
cirrhosis
CRC
CRC
CRC
CRC
CRC
CRC
CRC
CRC
CRC
CRC
hypertension
IBD
IBD
mucositis
peri-implantitis
schizophrenia
T1D
T1D
T2D
T2D
oralcavity
)
array6=(
rf
lsvm
svm
enet
lasso
)
array7=(
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
for index4 in ${!array6[*]}
do

if [ ${array6[index4]} == 'rf' ]
then

python classification.py ../data/shotgun/${array2[index2]}/${array3[index3]}/${array4[index1]}.txt ../results/classification/shotgun/${array7[index4]}/${array2[index2]}/${array3[index3]}/${array4[index1]}_${array6[index4]} -d 1:${array1[index1]}:${array5[index1]} -g [] -w -l ${array6[index4]} -nt 500 -nsl 5 -c entropy -nc 10

else

python classification.py ../data/shotgun/${array2[index2]}/${array3[index3]}/${array4[index1]}.txt ../results/classification/shotgun/${array7[index4]}/${array2[index2]}/${array3[index3]}/${array4[index1]}_${array6[index4]} -d 1:${array1[index1]}:${array5[index1]} -l ${array6[index4]}

fi

done
done
done
done

array1=(
rf
lsvm
svm
enet
lasso
)
array2=(
Random_forest
lsvm
svm
enet
lasso
)
array3=(
abundance
bool0
)

array4=(
edd_singh
crc_baxter
crc_zeller
crc_zhao
crc_xiang
ibd_gevers_2014
ibd_huttenhower
ibd_alm
ibd_engstrand_maxee
hiv_noguerajulian
hiv_dinh
hiv_lozupone
asd_son
autism_kb
t1d_alkanani
t1d_mejialeon
nash_chan
nash_baker
ra_littman
mhe_zhang
par_scheperjans
cdi_schubert
non_cdi_schubert
cdi_vincent_v3v5
cdi_youngster
ob_goodrich
ob_gordon_2008_v2
ob_zupancic
ob_ross
ob_baker
)
array5=(
EDD
CRC
CRC
CRC
CRC
CD
UC:CD
UC:CD
UC:CD
HIV
HIV
HIV
ASD
ASD
T1D
T1D
NASH
NASH
PSA:RA
CIRR:MHE
PAR
CDI
nonCDI
CDI
CDI
OB
OB
OB
OB
nonNASH-OB
)
array6=(
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
nash_wong
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
for index2 in ${!array3[*]}
do
for index3 in ${!array4[*]}
do

if [ ${array1[index1]} == 'rf' ]

then

python classification_16s.py ../data/16s/${array3[index2]}/${array4[index3]}.csv ../results/classification/16s/${array2[index1]}/${array3[index2]}/${array6[index3]}_${array1[index1]} -d 1:study_condition:${array5[index3]} -g [] -w -l ${array1[index1]} -nt 500 -nsl 5 -c entropy -nc 10

else

python classification_16s.py ../data/16s/${array3[index2]}/${array4[index3]}.csv ../results/classification/16s/${array2[index1]}/${array3[index2]}/${array6[index3]}_${array1[index1]} -d 1:study_condition:${array5[index3]} -l ${array1[index1]} 

fi

done
done
done

array1=(
rf
lsvm
svm
enet
lasso
)

array2=(
Random_forest
lsvm
svm
enet
lasso
)

array3=(
abundance
bool0
bool01
bool001
bool0001
bool00001
)

array4=(
species
genus
order
family
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

for index1 in ${!array1[*]}
do
for index2 in ${!array3[*]}
do
for index3 in ${!array4[*]}
do
for index4 in ${!array5[*]}
do

if [ ${array1[index1]} == 'rf' ]

then

python classification.py ../data/shotgun/lodo/${array3[index2]}/${array4[index3]}/profiles_metadata_table_crc.txt ../results/classification/shotgun/lodo/${array2[index1]}/${array3[index2]}/${array5[index4]}_crc_${array1[index1]} -d 1:study_condition:CRC -g [] -w -t dataset_name:${array5[index4]} -l ${array1[index1]} -nt 500 -nsl 5 -c entropy -nc 10

else

python classification.py ../data/shotgun/lodo/${array3[index2]}/${array4[index3]}/profiles_metadata_table_crc.txt ../results/classification/shotgun/lodo/${array2[index1]}/${array3[index2]}/${array5[index4]}_crc_${array1[index1]} -d 1:study_condition:CRC -t dataset_name:${array5[index4]} -l ${array1[index1]}

fi

done
done
done
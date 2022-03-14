array1=(
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

array2=(
species
genus
order
family
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
Schizophrenia
T1D
T1D
T2D
T2D
oralcavity
)

array5=(
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
control
stool
)

array6=(
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

for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do
for index3 in ${!array3[*]}
do


python species_test.py -d ../data/shotgun/${array3[index3]}/${array2[index2]}/${array1[index1]}.txt -o ../results/biomarkers/shotgun/${array3[index3]}/${array2[index2]}/${array1[index1]}.csv -c1 ${array4[index1]} -c2 ${array5[index1]} -t ${array6[index1]} -lvl ${array3[index3]}

done
done
done

array1=(
edd_singh
cdi_schubert
non_cdi_schubert
cdi_vincent_v3v5
cdi_youngster
ob_goodrich
ob_gordon_2008_v2
ob_zupancic
ob_ross
nash_ob_baker
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
nash_wong
nash_ob_baker
ra_littman
mhe_zhang
par_scheperjans
)

array2=(
edd_singh
cdi_schubert
non_cdi_schubert
cdi_vincent
cdi_youngster
ob_goodrich
ob_turnbaug
ob_zupancic
ob_ross
ob_zhu
crc_baxter
crc_zeller
crc_Wang
crc_chen
ibd_gevers_2014
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
)

array3=(
EDD
CDI
nonCDI
CDI
CDI
OB
OB
OB
OB
nonNASH-OB
CRC
CRC
CRC
CRC
CD
UC
UC
UC
HIV
HIV
HIV
ASD
ASD
T1D
T1D
NASH
NASH
PSA
CIRR
PAR
)

array4=(
abundance
bool0
)
array5=(
H
H
H
H
H
H
H
H
H
H
H
H
H
H
nonIBD
H
nonIBD
H
H
H
H
H
H
H
H
H
H
H
H
H
)



for index1 in ${!array1[*]}
do
for index2 in ${!array4[*]}
do

python species_test.py -d ../data/16s/${array4[index2]}/${array1[index1]}.csv -o ../results/biomarkers/16s/${array4[index2]}/${array2[index1]}.csv -c1 ${array3[index1]} -c2 ${array5[index1]} -t study_condition -lvl ${array4[index2]}

done
done




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

for index1 in ${!array1[*]}
do
for index2 in ${!array2[*]}
do

python species_test.py -d ../data/shotgun/lodo/${array1[index1]}/${array2[index2]}/profiles_metadata_table_crc.txt -o ../results/biomarkers/shotgun/lodo/${array1[index1]}/${array2[index2]}/species_crc.txt  -c1 CRC -c2 control -t study_condition -lvl ${array1[index1]}
done
done


#mkdir raw
#mv raw_* raw

#mkdir annot
#mv CoV2_genome* annot

#mkdir contour
#mv CoV2_contour* contour

mkdir ct_all
mkdir ct_B117
mkdir ct_P1
mkdir ct_B1351
mkdir ct_B1617
mkdir ct_Omicron

mv ct*_all*png ct_all
mv ct*_B117*png ct_B117
mv ct*_P1*png ct_P1
mv ct*_B1351*png ct_B1351
mv ct*_B1617*png ct_B1617
mv ct*_Omicron*png ct_Omicron

mkdir vgm
mv Vgm_* vgm

mkdir params
mv params_* params

mkdir krige_input
mkdir krige_output
mv kriged_* krige_output
mv krige_t* krige_input

mkdir hpc_logs
#mv *.R.o* hpc_logs/
mv *Rout hpc_logs/
mv *err.log hpc_logs/
#mv *.o* hpc_logs/

#mkdir minVar
#mv *.csv minVar

#tar -cvzf raw.tar.gz raw
#tar -cvzf annot.tar.gz annot

tar -cvzf ct_all.tar.gz ct_all
tar -cvzf ct_B117.tar.gz ct_B117
tar -cvzf ct_P1.tar.gz ct_P1 
tar -cvzf ct_B1351.tar.gz ct_B1351 
tar -cvzf ct_B1617.tar.gz ct_B1617
tar -cvzf ct_Omicron.tar.gz ct_Omicron

tar -cvzf krigeIn.tar.gz krige_input
tar -cvzf krigeOut.tar.gz krige_output
tar -cvzf params.tar.gz params

#tar -cvzf nsp7.tar.gz nsp7
#tar -cvzf nsp8.tar.gz nsp8
#tar -cvzf Pol.tar.gz Pol
#tar -cvzf Spike.tar.gz Spike

#tar -cvzf Pol_minVar.tar.gz Pol_minVar
#tar -cvzf Spike_minVar_maps.tar.gz Spike_minVar_maps

# cc plots for VOCs
cd /gpfs/home/loguerci/coronavirus/Covid19/Code/cc_strains

mv cc_B117*png B117
mv cc_P1*png P1
mv cc_B1351*png B1351
mv cc_B1617*png B1617

mkdir VOC_cc_lineplots
mv *png VOC_cc_lineplots

tar -cvzf cc_B117.tar.gz B117
tar -cvzf cc_P1.tar.gz P1 
tar -cvzf cc_B1351.tar.gz B1351 
tar -cvzf cc_B1617.tar.gz B1617

tar -cvzf VOC_cc_lineplots.tar.gz VOC_cc_lineplots









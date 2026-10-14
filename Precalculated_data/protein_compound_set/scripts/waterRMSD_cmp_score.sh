
for tg in $(cat ./pdbbind_exclude30.txt) ;
do
    echo $tg
	
    #python waterRMSD_score.py foldx ref/$tg_water.pdb FoldX/$tg.pdb > result_score/FoldX/$tg.dat
    python waterRMSD_score.py rism ref/${tg}_water.pdb 3D-RISM/$tg.pdb > result_score/3drism/$tg.dat
    #python waterRMSD_score.py wkgb ref/$tg_water.pdb wkgb/$tg.pdb > result_score/GalaxyWater-wKGB/$tg.dat
    python waterRMSD_score.py cnn ref/${tg}_water.pdb GalaxyWater-CNN/$tg.pdb > result_score/GalaxyWater-CNN/$tg.dat
    python waterRMSD_score.py gnn ref/${tg}_water.pdb WatGNN/${tg}_wat_lig_all.pdb > result_score/WatGNN/$tg.dat
    python waterRMSD_score.py hydraprot ref/${tg}_water.pdb HydraProt/${tg}_complex_waters.pdb > result_score/HydraProt/$tg.dat
    python waterRMSD_score.py WarPP ref/${tg}_water.pdb WarPP/${tg}_pred.pdb > result_score/WarPP/$tg.dat

done
#./summary_score.py pdbbind_exclude30.txt result_score/FoldX > summary_score/native_summary_FoldX.dat
./summary_score.py pdbbind_exclude30.txt result_score/3drism > summary_score/native_summary_3drism.dat
#./summary_score.py pdbbind_exclude30.txt result_score/GalaxyWater-wKGB > summary_score/native_summary_GalaxyWater-wKGB.dat
./summary_score.py pdbbind_exclude30.txt result_score/GalaxyWater-CNN > summary_score/native_summary_GalaxyWater-CNN.dat
./summary_score.py pdbbind_exclude30.txt result_score/WatGNN > summary_score/native_summary_WatGNN.dat
./summary_score.py pdbbind_exclude30.txt result_score/HydraProt > summary_score/native_summary_HydraProt.dat
./summary_score.py pdbbind_exclude30.txt result_score/WarPP > summary_score/native_summary_WarPP.dat
python ./plot_pdbbind_score.py

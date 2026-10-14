
for tg in $(cat ./mobywat_trg.txt) ;
do
    echo $tg
	
    python  waterRMSD_score.py rism      water_35/$tg.pdb water_cutoff_3.5/rism_mobywat/${tg}_pred.pdb      > result_35_score/3drism/$tg.dat
    #python waterRMSD_score.py foldx     water_35/$tg.pdb water_cutoff_3.5/FoldX/${tg}_pred.pdb             > result_35_score/FoldX/$tg.dat
    python  waterRMSD_score.py cnn       water_35/$tg.pdb water_cutoff_3.5/CNN_mobywat/${tg}_pred.pdb       > result_35_score/GalaxyWater-CNN/$tg.dat
    #python waterRMSD_score.py wkgb      water_35/$tg.pdb water_cutoff_3.5/wkgb/${tg}_pred.pdb              > result_35_score/GalaxyWater-wKGB/$tg.dat
    python  waterRMSD_score.py gnn       water_35/$tg.pdb water_cutoff_3.5/GNN_mobywat/${tg}_pred.pdb       > result_35_score/WatGNN/$tg.dat
    python  waterRMSD_score.py hydraprot water_35/$tg.pdb water_cutoff_3.5/hydraprot_mobywat/${tg}_pred.pdb > result_35_score/HydraProt/$tg.dat
    python  waterRMSD_score_mobywat.py mobywat   water_35/$tg.pdb water_cutoff_3.5/mobywat_mtol15/${tg}_pred.pdb    > result_35_score/MobyWat/$tg.dat
    #python waterRMSD_score_superwater.py water_35/$tg.pdb water_cutoff_3.5/Superwater/superwater_010/${tg}_pred.pdb > result_35_score/SuperWater/$tg.dat
	
done


./summary_native_score.py mobywat_trg_5.txt  result_35_score/3drism           > summary_35_score/native_summary_3drism.dat
#./summary_native_score.py mobywat_trg_5.txt result_35_score/FoldX            > summary_35_score/native_summary_FoldX.dat
./summary_native_score.py mobywat_trg_5.txt  result_35_score/GalaxyWater-CNN  > summary_35_score/native_summary_GalaxyWater-CNN.dat
#./summary_native_score.py mobywat_trg_5.txt result_35_score/GalaxyWater-wKGB > summary_35_score/native_summary_GalaxyWater-wKGB.dat
./summary_native_score.py mobywat_trg_5.txt  result_35_score/WatGNN           > summary_35_score/native_summary_WatGNN.dat
./summary_native_score.py mobywat_trg_5.txt  result_35_score/HydraProt        > summary_35_score/native_summary_HydraProt.dat
#./summary_native_score.py mobywat_trg_5.txt result_35_score/SuperWater       > summary_35_score/native_summary_SuperWater.dat
./summary_native_score.py mobywat_trg_5.txt  result_35_score/MobyWat          > summary_35_score/native_summary_MobyWat.dat

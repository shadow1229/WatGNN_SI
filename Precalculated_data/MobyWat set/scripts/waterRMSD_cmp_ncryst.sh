
for tg in $(cat ./mobywat_trg.txt) ;
do
    echo $tg

    python  waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/rism_mobywat/${tg}_pred.pdb      > result_35/3drism/$tg.dat
    #python waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/FoldX/${tg}_pred.pdb             > result_35/FoldX/$tg.dat
    python  waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/CNN_mobywat/${tg}_pred.pdb       > result_35/GalaxyWater-CNN/$tg.dat
    #python waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/wkgb/${tg}_pred.pdb              > result_35/GalaxyWater-wKGB/$tg.dat
    python  waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/GNN_mobywat/${tg}_pred.pdb       > result_35/WatGNN/$tg.dat
    python  waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/hydraprot_mobywat/${tg}_pred.pdb > result_35/HydraProt/$tg.dat
    python  waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/mobywat_mtol15/${tg}_pred.pdb    > result_35/MobyWat/$tg.dat
    #python waterRMSD_ncryst.py water_35/$tg.pdb water_cutoff_3.5/Superwater/superwater_010/${tg}_pred.pdb > result_35/SuperWater/$tg.dat
	
done

./summary_native.py mobywat_trg_5.txt  result_35/3drism           > summary_35/native_summary_3drism.dat
#./summary_native.py mobywat_trg_5.txt result_35/FoldX            > summary_35/native_summary_FoldX.dat
./summary_native.py mobywat_trg_5.txt  result_35/GalaxyWater-CNN  > summary_35/native_summary_GalaxyWater-CNN.dat
#./summary_native.py mobywat_trg_5.txt result_35/GalaxyWater-wKGB > summary_35/native_summary_GalaxyWater-wKGB.dat
./summary_native.py mobywat_trg_5.txt  result_35/WatGNN           > summary_35/native_summary_WatGNN.dat
./summary_native.py mobywat_trg_5.txt  result_35/HydraProt        > summary_35/native_summary_HydraProt.dat
#./summary_native.py mobywat_trg_5.txt result_35/SuperWater       > summary_35/native_summary_SuperWater.dat
./summary_native.py mobywat_trg_5.txt  result_35/MobyWat          > summary_35/native_summary_MobyWat.dat

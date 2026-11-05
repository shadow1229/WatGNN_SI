
for tg in $(cat ./targets_cmp) ;
do
    echo $tg

    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb 3drism/$tg.pdb > result/3drism/$tg.dat
    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb FoldX/$tg.pdb > result/FoldX/$tg.dat
    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb GalaxyWater-CNN/$tg.pdb > result/GalaxyWater-CNN/$tg.dat
    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb wkgb/$tg.pdb > result/GalaxyWater-wKGB/$tg.dat
    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb watGNN/$tg.pdb > result/WatGNN/$tg.dat
    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb Hydraprot_sorted/$tg.pdb > result/HydraProt/$tg.dat
    #python waterRMSD_ncryst.py wkgb_ref/$tg.pdb Superwater/superwater_010/$tg.pdb > result/SuperWater/$tg.dat
    python waterRMSD_ncryst.py wkgb_ref/$tg.pdb MobyWat/MobyWat_aligned/$tg.pdb     > result/MobyWat/$tg.dat
    python waterRMSD_ncryst.py wkgb_ref/$tg.pdb MobyWat/MobyWat_aligned_1ns/$tg.pdb > result/MobyWat_1ns/$tg.dat
done

#./summary_native.py targets_cmp result/3drism > summary/native_summary_3drism.dat
#./summary_native.py targets_cmp result/FoldX > summary/native_summary_FoldX.dat
#./summary_native.py targets_cmp result/GalaxyWater-CNN > summary/native_summary_GalaxyWater-CNN.dat
#./summary_native.py targets_cmp result/GalaxyWater-wKGB > summary/native_summary_GalaxyWater-wKGB.dat
#./summary_native.py targets_cmp result/WatGNN > summary/native_summary_WatGNN.dat
#./summary_native.py targets_cmp result/HydraProt > summary/native_summary_HydraProt.dat
#./summary_native.py targets_cmp result/SuperWater > summary/native_summary_SuperWater.dat
./summary_native.py targets_cmp result/MobyWat     > summary/native_summary_MobyWat.dat
./summary_native.py targets_cmp result/MobyWat_1ns > summary/native_summary_MobyWat_1ns.dat

python ./plot_wkgb_ncryst.py
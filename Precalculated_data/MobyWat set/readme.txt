original_PDB     : folders for original PDB file for mobywat set
protein          : Exclude water from original_PDB, used as inputs for WatGNN, GalaxyWater_CNN,3D-RISM,HydraProt
                   (3D-RISM uses scripts/build_mobywat_input_script/prepare_placevent_inputs.py for setup procedure including adding hydrogen)
water_35 : water molecules from the crystal structure, where the water molecule has B-Factor<40 and 
                   having distance closer than 3.5 Angstrom from both receptor and ligand.
scripts          : scripts for building the input structure, performance evaluation, and plotting

GNN_mobywat      : WatGNN prediction folder ( saved as [PDB ID]_pred.pdb )
CNN_mobywat      : GalaxyWater-CNN prediction folder ( saved as [PDB ID].pdb )
rism_mobywat     : 3D-RISM prediction folder ( saved as [PDB ID]_pred.pdb )
hydraprot_mobywat: HydraProt prediction folder ( saved as [PDB ID]_waters.pdb )
mobywat_method3  : Mobyat prediction folder, from MobyWat's publication and the best method, method 3 was selected.
 - O_*_mtIDa.lst	match list
 - O_*_prIDa.pdb	prediction list
 - O_*_rIFw.pdb	    reference waters
water_cutoff_3.5 : water molecules from the prediction methods, where the predicted position 
                   having distance closer than 4.5 (3.5+1) Angstrom from both receptor and ligand.
result       : performance result for each target for fixed prediction budget
result_score : performance result for each target for fixed score cutoff
summary      : 5 target-averaged result for each target for fixed prediction budget
summary_score: 5 target-averaged result for each target for fixed score cutoff
plots        : plots used in the manuscript's Supplementary information